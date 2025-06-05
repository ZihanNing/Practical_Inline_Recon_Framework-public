function svr=svrAlternateMinimization(svr,nIt,tol)

%SVRALTERNATEMINIMIZATION   Performs an alternate minimization for SVR
%   SVR=SVRALTERNATEMINIMIZATION(SVR,{NIT},{TOL})
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   {NIT} is the maximum number of iterations
%   {TOL} is the tolerance set for convergence
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)

if ~exist('nIt','var') || isempty(nIt);nIt=300;end%300;end
if ~exist('tol','var') || isempty(tol);tol=1e-3;end%5e-3;end%1e-2;end

%HARDCODED MAXIMUM NUMBER OF ITERATIONS FOR EACH MOTION ESTIMATE STEP
maxItMot=[128 64 8]*4;
%maxItMot=[0 0 10];
nItMot=[0 0 0];
nItReset=[1 1 1];
nItResetTotal=[1 1 1];
%if svr.PreTr;maxItMot=maxItMot/2;end

xPrev=svr.x;
if isempty(svr.xx)
    %SOFT MASK
    svr=svrGenerateSoftMask(svr);

    %WRITE AVERAGED INFORMATION FOR DEBUGGING
    svr=svrEmpiricalPseudoInverseAveraging(svr);

    %SEGMENTATION
    svr=svrRunSegmentation(svr,0);
end

%GENERATE THE MASK USED IN REGISTRATION
svr=svrGenerateSoftMask(svr,2);

%PARAMETERS FOR ROBUST RECONSTRUCTION
redN=0.5;%Reduction factor to approach the target norm
%redN=0.925;%Reduction factor to slowly approach the target norm

conv=1;
if ~svr.PreTr;motEst=0;elseif svr.PreTr==1;motEst=max(svr.ParSVR.EstT-1,0);else motEst=max(svr.ParSVR.EstT+1,0);end
n=1;
if ~svr.PreTr;svr.x(:)=0;else svr.x=xPrev;end
xPrev=[];

while(1)    
    %if conv==1 && motEst~=svr.ParSVR.EstT;non=2;end%Initial homotopy value
    %if motEst==svr.ParSVR.EstT;non=svr.ParSVR.Lp;end
    if ~svr.PreTr;non=2;else non=svr.ParSVR.Lp;end
    
    %%WE COMPUTE THE TRACES    
    svr.xx=plugNoise(svr.x,1);    
    svr.xx=bsxfun(@rdivide,svr.xx,sqrt(normm(svr.xx)/numel(svr.xx)));
    svr=svrEncode(svr,0);
    trFid=0;
    for v=1:svr.NV
        if ~isfield(svr,'EnViews') || svr.EnViews(v)
            %SLICE WEIGHTS
            svr.yy{v}=bsxfun(@times,svr.yy{v},sqrt(svr.W{v}));
            %APODIZATION
            if svr.ParSVR.UseApo>0;svr.yy{v}=bsxfun(@times,svr.yy{v},sqrt(svr.H{v}));end
            trFid=trFid+normm(svr.yy{v});            
        end
    end
    
    %WE COMPUTE THE REGULARIZATION WEIGHTS
    if svr.ParSVR.tiSh~=0
        svr.We=abs(shearletTransformGPUMemory(svr.x,svr.Sh));
        svr.We=1./(svr.We+1e-9);
        svr.xx=plugNoise(svr.x,1);
        svr.xx=bsxfun(@rdivide,svr.xx,sqrt(normm(svr.xx)/numel(svr.xx)));
        svr.xx=shearletTransformGPUMemory(svr.xx,svr.Sh);
        svr.xx=svr.xx.*sqrt(svr.We);
        trReg=normm(svr.xx);
        svr.ParSVR.tiSh=gather(20*2*trFid/(trReg+eps));
        fprintf('Regularization factor: %.6g\n',svr.ParSVR.tiSh);
    end
    
    %visReconstruction(svr.x)
    svr=svrCG(svr,1+(n<2)+(n<3),[],[],7);
    %visReconstruction(svr.x)
    
    

    %if motEst==2
    %    visReconstruction(svr.x)
    %end
    
    non=max(svr.ParSVR.Lp,non-redN); 
    if svr.ParSVR.Debug
        if ~isfield(svr,'xHist');svr.xHist=gather(svr.x);else svr.xHist=cat(4,svr.xHist,gather(svr.x));end
    end
        
    if conv==1 && motEst<=svr.ParSVR.EstT && (motEst<svr.ParSVR.EstT || (non==svr.ParSVR.Lp && motEst==svr.ParSVR.EstT))
        svr=svrCNNFinalBrainDetection(svr);        
        %if motEst==0;svr=svrGenerateSoftMask(svr,4);end%TO USE THE CNN MASK, NOT ROBUST ENOUGH
        svrWriteData(svr,sprintf('Step%02d-After%dMot-CNNMask',svr.ParSVR.Step,motEst),svr.Elx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;        
        %svrWriteData(svr,sprintf('Step%02d-After%dMot-RecMask',svr.ParSVR.Step,motEst),svr.M,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;        
        svrWriteData(svr,sprintf('Step%02d-After%dMot',svr.ParSVR.Step,motEst),svr.x,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;
        svr.yy=svr.y;
        svr=svrDecode(svr,1);
        xx=svr.xx;
        if motEst~=0
            svrWriteData(svr,sprintf('Step%02d-StacksSpatialCoordinates',svr.ParSVR.Step),xx,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;           
            svrWriteData(svr,sprintf('Step%02d-SpatialCoordinates',svr.ParSVR.Step),mean(xx,4),[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;   
        end
    end

    svr=svrSliceWeights(svr,non);       
    
    %if (motEst==svr.ParSVR.EstT && non==svr.ParSVR.Lp) || n>nIt;break;end     
    if motEst==svr.ParSVR.EstT || n>nIt;break;end     

    if motEst==0 && svr.ParSVR.EstT>0
        if nItMot(1)<maxItMot(1);[svr,conv]=svrSolveTVolu(svr);
        else conv=1;svr.w=[];svr.convT=[];svr=rmfield(svr,{'w','convT'});
        end        
    elseif motEst==1
        if nItMot(2)<maxItMot(2);[svr,conv]=svrSolveTPack(svr);            
        else conv=1;svr.w=[];svr.convT=[];svr=rmfield(svr,{'w','convT'});
        end
    elseif motEst==2
        if nItMot(3)<maxItMot(3);[svr,conv]=svrSolveTExci(svr);
        else conv=1;svr.w=[];svr.convT=[];svr=rmfield(svr,{'w','convT'});
        end        
    end 
    if motEst<=2
        m=motEst+1;
        nItMot(m)=nItMot(m)+1;motEst=motEst+conv;
        if nItMot(m)==nItResetTotal(m)
            if isfield(svr,'convT');svr=rmfield(svr,'convT');end
            nItReset(m)=nItReset(m)+1;
            nItResetTotal(m)=nItResetTotal(m)+nItReset(m);            
        end
    end
    n=n+1;
end
%toc
if n==nIt;fprintf('Alternate minimization solver terminated without reaching convergence\n');end

%WE WRITE THE DATA
if svr.PreTr;svrWriteData(svr,sprintf('Step%02d-IntermediateResult',svr.ParSVR.Step),svr.x,[],0.8,[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;end%WE WRITE THE DATA CROPPING TO POSITIVE AND AT 0.8

if svr.ParSVR.Debug
    %WE WRITE THE CONVERGENCE HISTORY
    svrWriteData(svr,sprintf('History-Reconstructions%.2f',1/mean(svr.MSS)),svr.xHist,[],[],[],1);
    svr=rmfield(svr,'xHist');

    %WE WRITE THE RESIDUALS
    svrWriteData(svr,sprintf('History-Residuals%.2f',1/mean(svr.MSS)),svr.E,1,[],[],1);
end
