function svr=svrCNNBrainDetection(svr)

%SVRCNNBRAINDETECTION   Detects the brain based on V-Net CNN training
%   SVR=SVRCNNBRAINDETECTION(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%

gpu=svr.rec{1}.Dyn.GPU;
svr.NV=length(svr.y);
cnnPath=strcat(fileparts(mfilename('fullpath')),'/../../../../Software/CNNSVR');
%cnnPath='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/CNNSVR';
networkDir=sprintf('%s/Network',cnnPath);%Placeholder for training images

load(sprintf('%s/fetalbrainnettcGross%02d-%d.mat',networkDir,1,3),'fetalbrainnettc');
load(sprintf('%s/fetalbrainnettc%02d-%d-Res4.mat',networkDir,1,3),'fetalbrainnettd');

networkRes=1.2;%Resolution of the data used to train the network
NChan=4;
NN=[16 10 8]*16;%Standard volume size
Res=4;%Resolution to obtain the segmentations
NNr=NN/Res;%Standard volume size at the resolution to obtain the segmentations
NNr2=[NNr(1) prod(NNr(2:3)) NChan];%Size of the input to the network
skipDeconv=1;

NChanCrop=4;
NNCrop=[8 8 8]*16;%Standard volume size
ResCrop=4;%Resolution to obtain the segmentations
NNrCrop=NNCrop/ResCrop;%Standard volume size at the resolution to obtain the segmentations
NNr2Crop=[NNrCrop(1) prod(NNrCrop(2:3)) NChanCrop];%Size of the input to the network
skipDeconvCrop=1;

for v=1:svr.NV    
    %LOWRES ESTIMATE
    x=svr.y{v};
    if gpu;x=gpuArray(x);end
    NX=size(x);
    NXRes=round(NX.*svr.MS{v}/networkRes);
    x=resampling(x,NXRes,[],2*ones(1,3));      
    x=resampling(x,NN,2);
    xM=multDimMea(abs(x),[1:2 4]);
    if NChan==8;multZ=2;else multZ=1;end
    xM=reshape(xM,[Res/multZ multZ*NN(3)/Res]);
    [~,iMax]=max(xM,[],1);
    nMax=(1:Res/multZ:NN(3))+iMax-1;        
    x=x(:,:,nMax,:);
    if NChan~=1;x=resampling(x,NNr(1:2)*2,0,ones(1,2));end
    if NChan==4
        x=reshape(x,[2 NNr(1) 2 NNr(2) NNr(3)]);
        x=permute(x,[2 4 5 1 3]);
        x=reshape(x,[NNr 4]);
    elseif NChan==8
        x=reshape(x,[2 NNr(1) 2 NNr(2) 2 NNr(3)]);
        x=permute(x,[2 4 6 1 3 5]);
        x=reshape(x,[NNr 8]);        
    else
        x=resampling(x,NNr,0,ones(1,3));
    end
    x=reshape(x,NNr2);
    x=gather(single(abs(x)));
    y=predict(fetalbrainnettc,x)';
    y=reshape(y,NNr/2^skipDeconv);
    y=resampling(y,NNr,[],2*ones(1,3));
    zlr=(centerMass(max(y,0).^2)-0.5)*Res+0.5;
    
    %HIGHRES ESTIMATE
    x=svr.y{v};
    if gpu;x=gpuArray(x);end
    NX=size(x);
    NXRes=round(NX.*svr.MS{v}/networkRes);
    x=resampling(x,NXRes,[],2*ones(1,3));
    x=resampling(x,NN,2);
    zr=round(zlr);
    H=cell(1,3);
    zrc=zr-NN/2;
    for n=1:3;H{n}=-zrc(n);end
    x=shifting(x,H);
    x=resampling(x,NNCrop,2);    
    
    xM=multDimMea(abs(x),[1:2 4]);
    if NChan==8;multZ=2;else multZ=1;end
    xM=reshape(xM,[ResCrop/multZ multZ*NNCrop(3)/ResCrop]);
    [~,iMax]=max(xM,[],1);
    nMax=(1:ResCrop/multZ:NNCrop(3))+iMax-1;        
    x=x(:,:,nMax,:);
    if NChanCrop~=1;x=resampling(x,NNrCrop(1:2)*2,0,ones(1,2));end
    if NChanCrop==4
        x=reshape(x,[2 NNrCrop(1) 2 NNrCrop(2) NNrCrop(3)]);
        x=permute(x,[2 4 5 1 3]);
        x=reshape(x,[NNrCrop 4]);
    elseif NChanCrop==8
        x=reshape(x,[2 NNrCrop(1) 2 NNrCrop(2) 2 NNrCrop(3)]);
        x=permute(x,[2 4 6 1 3 5]);
        x=reshape(x,[NNrCrop 8]);        
    else
        x=resampling(x,NNrCrop,0,ones(1,3));
    end
    x=reshape(x,NNr2Crop);
    x=gather(single(abs(x)));
    y=predict(fetalbrainnettd,x)';
    y=reshape(y,NNrCrop/2^skipDeconv);
    y=resampling(y,NNrCrop,[],2*ones(1,3));    
    
    zhr=(centerMass(max(y,0).^4)-0.5)*Res+0.5;
    zhr=zhr+(NN-NNCrop)/2+zrc+ceil((NXRes-NN)/2);   
    %RADII
    if gpu;y=gpuArray(y);end
    y=gather(resampling(y,NNCrop,[],2*ones(1,3)));
    y=y>0.5;
    if any(y(:)==1)
        ROI=computeROI(gather(y>0.5));       
        zhr(4:6)=diff(ROI(1:3,1:2),1,2)/2;    
    else
        zhr(4:6)=zeros(1,3);
    end
    w=buildGeometry(NXRes,'Ellipsoid',[zhr zeros(1,3)],gpu);
    dist=4;%THIS IS A HARDCODED PARAMETER       
    w=morphFourier(w,dist*ones(1,3),svr.MS{v},ones(1,3),1);
    w=resampling(w,NX,[],2*ones(1,3));
    w=single(w>0.5);
    svr.El{v}=gather(w);
    zhr=[(zhr(1:3)-0.5).*NX./NXRes+0.5 zhr(4:6).*NX./NXRes];    
    svr.Elp{v}=zhr;
    svr.y{v}=svr.y{v};
end
if svr.ParSVR.Visualize;svrVisualization(svr.y,svr.El,'CNN brain detection');end
svrWriteData(svr,'El',svr.El,1,[],'',0);
