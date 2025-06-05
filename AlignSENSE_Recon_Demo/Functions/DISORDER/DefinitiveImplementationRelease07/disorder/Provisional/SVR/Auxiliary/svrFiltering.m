function svr=svrFiltering(svr,typ)

%SVRFILTERING   Performs some filtering of the original images
%   SVR=SVRSETUP(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2;typ=1;end

gpu=svr.rec{1}.Dyn.GPU;
svr.NV=length(svr.y);

svr.yOr=svr.y;%Unfiltered data
if typ==1
    if svr.ParSVR.GibbsRingi~=0
        for v=1:svr.NV
            svr.yOr{v}=gather(svr.yOr{v});
            NY=size(svr.y{v});NY(end+1:3)=1;NY=NY(1:2);
            H=buildFilter(NY,'tukeyIso',ones(1,2),gpu,svr.ParSVR.GibbsRingi);
            y=svr.y{v};
            if gpu;y=gpuArray(y);end        
            y=filtering(y,H);
            svr.y{v}=gather(y);
        end
    end
else
    J=1;
    bhat=1;%0.5;%0.5;%[];%0.05;
    %wi=[16 32];
    wi=[8 32];
    for v=1:svr.NV
        svr.yOr{v}=gather(svr.yOr{v});
        NY=size(svr.y{v});NY(end+1:3)=1;
        sH=buildShearlet(NY(1:2),J,gpu,ones(1,J),'kos');
        x=svr.y{v};
        if isfield(svr,'no');nF=svr.no{v};
        else nF=real(x);nF(:)=1;
        end
        if gpu;[x,nF]=parUnaFun({x,nF},@gpuArray);end      
        x=shearletFilter(x,sH,nF,bhat,wi);
        svr.y{v}=gather(x);
        if isfield(svr,'no');svr.no{v}=gather(svr.no{v});end
    end   
    %WE WRITE THE FILTERED DATA
    svrWriteData(svr,'Re',svr.y,1,[],'',0);
end