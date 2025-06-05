function svr=svrCorrectInhom(svr,dynR,filtS)

%SVRCORRECTINHOM   Corrects inhomogeneity of images based on a B1 field
%   SVR=SVRCORRECTINHOM(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), 
%   orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2 || isemtpy(dynR);dynR=0.1;end
if nargin<3 || isemtpy(filtS);filtS=0.1;end%0.05 before

gpu=svr.rec{1}.Dyn.GPU;

svr.NV=length(svr.y);

%MAP THE COORDINATES BACK
for v=1:svr.NV
    %WE COMPUTE THE B1 FIELD
    %M=single(dynInd(svr.yB1{v},2,4)>1e-6);%For if we would like to use a mask
    xFA=abs(svr.yB1{v});%ABS here is important!
    y=svr.y{v};
    no=[];if isfield(svr,'no');no=svr.no{v};end
    if gpu;[xFA,y,no]=parUnaFun({xFA,y,no},@gpuArray);end  
    xFA=convertRotation(atan(2*dynInd(xFA,2,4)./(dynInd(xFA,1,4)+eps)),'rad','deg');%FA from DREAM paper
    xFA=(xFA./svr.recB1{v}.Par.Labels.FlipAngle(2)).^svr.ParSVR.CorrectInhomPowe;%FA normalized
    %svr.yB1{v}=gather(svr.yB1{v});

    %MAP VOLUME
    xFANY=mapVolume(xFA,y,svr.MTB1{v},svr.MT{v});
        
    %FILTERING
    N=size(xFANY);N=N(1:3);
    H=buildFilter(2*N,'tukeyIso',filtS*ones(1,3),gpu,1,1);    
    xFANY=filtering(xFANY,H,1);
    xFANY=max(xFANY,dynR);
    %svrVisualization(y,[],'Original data')
    y=y./abs(xFANY);
    if isfield(svr,'no');no=no./abs(xFANY);end
    %svrVisualization(y,[],'Corrected data')
    svr.y{v}=gather(y);
    if isfield(svr,'no');svr.no{v}=gather(no);end
    svr.yB1{v}=[];svr.recB1{v}=[];
end
svr=rmfield(svr,{'yB1','MSB1','MTB1','recB1'});
if svr.ParSVR.Visualize;svrVisualization(svr.y,[],'Inhomogeneity corrected');end

%WE WRITE THE INHOMOGENEITY CORRECTED DATA
svrWriteData(svr,'Ic',svr.y,1,[],'',0);
