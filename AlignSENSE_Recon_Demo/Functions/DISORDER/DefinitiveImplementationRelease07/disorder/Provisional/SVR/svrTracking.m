function svr=svrTracking(svr,typ)

%SVRTRACKING   Detects and tracks the brain on a per-volume basis
%   SVR=SVRTRACKING(SVR)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

gpu=isa(svr.x,'gpuArray');

%SOFT MASK
svr=svrGenerateSoftMask(svr);

%WRITE AVERAGED INFORMATION FOR DEBUGGING
svr=svrEmpiricalPseudoInverseAveraging(svr);

%SEGMENTATION
svr=svrRunSegmentation(svr,0);


if typ==0%TRANSLATIONAL TRACKING
    
    %%COMMENTING THIS PERFORMS THE REGISTRATION USING THE ELLIPSOIDS!
    svr.yy=svr.y;
    svr=svrDecode(svr,1,1);  
    
    %GENERATE THE MASK USED IN TRACKING
    %svr=svrGenerateSoftMask(svr,3);
    svr.M=svr.Mx;
    
    %TO ACCELERATE TRACKING
    dilate=60;
    svr.Mtr=morphFourier(single(svr.M>0),dilate,svr.MSS,ones(1,3));

    %ROI COMPUTATION
    ROI=computeROI(svr.Mtr);
    fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',ROI'));

    %ROI EXTRACTION FOR POSE ESTIMATION
    xx=extractROI(svr.xx,ROI,1,1:3);
    Mx=extractROI(svr.M,ROI,1,1:3);

    %REGISTRATION
    Tr=single(zeros([ones(1,3) svr.NV 6]));
    %kmax=[8 7 6 5];
    %dk=[4 3 2 1];
    %phmax=[40 40 40 40];
    %Nori=[0 0 0 0];
    %Nrot=[0 0 0 0];    
    %lev=[5 4 3 2];
    
%     kmax=[6 5];
%     dk=[2 1];
%     phmax=[40 40];
%     Nori=[0 0];
%     Nrot=[0 0];    
%     lev=[3 2];
    
    kmax=5;
    dk=1;
    phmax=40;
    Nori=0;
    Nrot=0;    
    lev=2;

    %svrVisualization(num2cell(xx,1:3),[],sprintf('Entry to translational tracking'),round(mean(svr.MSS(:))));

    
    BlSz=[100 10];
    metric='NC';
    fracOrder=0;
    useBest=1;
    [~,Tr,~]=integerShiftRotationRegistration(xx,Mx,Tr,kmax,dk,phmax,Nrot,Nori,lev,fracOrder,svr.MSS,BlSz,useBest,metric);xx=[];

    Tr=permute(Tr,[1 5 4 2 3]);
    Tr=Tr(1,1:3,:);%Only translations
    Tr=bsxfun(@minus,Tr,mean(Tr,3));%Demeaning to keep the center of the FOV
    Tr=bsxfun(@rdivide,Tr,svr.MSS);%We convert to pixels
    Tr=bsxfun(@plus,Tr,svr.cFOVRec-1);%We add to the center of the FOV
    Tr=permute(Tr,[2 1 3]);
    Tr(4,:,:)=1;
    Tr=matfun(@mldivide,svr.AcqToRec,Tr);%We obtain the new centers in the acquisition space
    Tr=Tr(1:3,:,:)+1;
    Tr=permute(Tr,[2 1 3]);
    for v=1:svr.NV
        fprintf('CNN brain centers stack %d:%s\n',v,sprintf(' %.2f',svr.Elp(1,1:3,v)));
        fprintf('Tracking brain centers stack %d:%s\n',v,sprintf(' %.2f',Tr(1,:,v)));
    end

    %WE ASSIGN THE NEW CENTERS
    svr.Elp(1,1:3,:)=Tr;
    
else%ROTATIONAL TRACKING
    svr.yy=svr.y;    
    svr=svrDecode(svr,2,1,1);%FOR HIGH INTENSITY SLICE SELECTION
    %svr=svrDecode(svr,2,1,0);

    
    %GENERATE THE MASK USED IN TRACKING
    %svr=svrGenerateSoftMask(svr,3);
    svr.M=svr.Mx;
    
%    lev=[5 4 3 2];
%    kmax=[2 2 2 2];
%    dk=[4 3 2 1];
%    %phmax=[60 45 30 15];       
%    %phmax=[90 70 50 30]*2;
%    phmax=[180 150 120 90]/3;%/3;
%    Nori=[30 25 20 15]*2;%Nori=[180 160 140 120];%In case more are needed
%    Nrot=[40 30 20 10]*2;%Nrot=[60 50 40 30];%In case more are needed    
    
    lev=[4 3 2];
    kmax=[2 2 2];
    dk=[3 2 1];
    %phmax=[60 45 30 15];       
    %phmax=[90 70 50 30]*2;
    phmax=[150 120 90]/3;%/3;
    Nori=[25 20 15]*2;%Nori=[180 160 140 120];%In case more are needed
    Nrot=[30 20 10]*2;%Nrot=[60 50 40 30];%In case more are needed    
    
    BlSz=[50 10];
    useBest=3;
    fracOrder=0;
    metric='NC';
    
    %TO ACCELERATE POSE ESTIMATION
    dilate=30;
    svr.Mtr=morphFourier(single(svr.M>0),dilate,svr.MSS,ones(1,3));

    %ROI COMPUTATION
    ROI=computeROI(svr.Mtr);
    fprintf('ROI processing:\n%s',sprintf(' %d %d %d %d %d %d\n',ROI'));
    
    %ROI EXTRACTION FOR POSE ESTIMATION
    xx=extractROI(svr.xx,ROI,1,1:3);
    Mx=extractROI(svr.M,ROI,1,1:3);
    center=svr.cFOVRec-ROI(:,5)';
    
    [~,svr.TV]=integerShiftRotationRegistration(xx,Mx,svr.TV,kmax,dk,phmax,Nrot,Nori,lev,fracOrder,svr.MSS,BlSz,useBest,metric,center);           
      
    %WE DETECT OUTLIERED VIEWS
    svr.yy=svr.My;
    svr=svrDecode(svr,1);   
    MMxx=mapToZeroOne(svr.xx,0.05);
    
    svr.EnViewsOr=multDimSum(svr.Mx.*MMxx,1:3);
    svr.EnViews=svr.EnViewsOr(:)';
    if svr.NV>svr.ParSVR.maxNoViews        
        [~,iS]=sort(svr.EnViews);
        svr.EnViews(iS(1:svr.NV-svr.ParSVR.maxNoViews))=0;
    end
    svr.EnViews=svr.EnViews/median(svr.EnViews);
    svr.EnViews=svr.EnViews>0.5;%LESS THAN HALF OVERLAP
    indDisc=find(~svr.EnViews);
    if ~isempty(indDisc);fprintf('DISCARDING STACKS%s FOR RECONSTRUCTION\n',sprintf(' %d',indDisc));end
        
    %SOFT MASK
    svr=svrGenerateSoftMask(svr);

    %WRITE AVERAGED INFORMATION FOR DEBUGGING
    svr=svrEmpiricalPseudoInverseAveraging(svr);

    %SEGMENTATION
    svr=svrRunSegmentation(svr,0);
end


