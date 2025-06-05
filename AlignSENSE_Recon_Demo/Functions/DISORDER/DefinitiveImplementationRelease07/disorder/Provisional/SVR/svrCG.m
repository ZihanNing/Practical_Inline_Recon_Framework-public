function svr=svrCG(svr,nIt,tolType,tol,NitTest,norm,initStep)

%SVRCG   Performs a CG-based pseudoinverse reconstruction for SVR
%   SVR=SVRCG(SVR,{NIT},{TOLTYPE},{TOL})
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   {NIT} is the maximum number of iterations
%   {TOLTYPE} is the type of tolerance used for convergence
%   {TOL} is the tolerance set for convergence
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)

if nargin<2 || isempty(nIt);nIt=300;end%300;end
if nargin<3 || isempty(tolType);tolType='Energy';end%tolType='NormwiseBackward2Error'
if nargin<4 || isempty(tol);tol=1e-2;end%5e-3;end%1e-2;end
if nargin<5 || isempty(tol);NitTest=1;end%5e-3;end%1e-2;end
if nargin<6 || isempty(norm);norm=0;end%5e-3;end%1e-2;end

gpu=isa(svr.x,'gpuArray');
tolTyp=stopCondition(tolType,'CG');

if svr.ParSVR.usePrecond==1%It makes no difference
    svr.xx=svr.x;svr.xx(:)=1;
    svr=svrDecode(svrEncode(svr));
    svr.Precond=(svr.xx+median(svr.xx(:))).^(-1);
end

svr.yy=svr.y;
if norm==1
    for v=1:svr.NV;svr.yy{v}=1;end
end
svr=svrDecode(svr);
r=svr.xx;

if nargin<7;initStep=all(svr.x(:)==0);end
if norm==1
    if initStep;svr.x(:)=0;else svr.x(:)=1;end
end
%initStep
svr.xx=svr.x;
%if norm==2;svr.xx(:)=0;end
svr=svrDecode(svrEncode(svr));
EHE=svr.xx;
for g=1:length(svr.ParSVR.regFracOrd)
    if svr.ParSVR.ti(g)~=0;EHE=EHE+svr.ParSVR.ti(g)*real(filtering(svr.x,svr.F{g},1));end%%%NOTE THAT THIS SHOULD BE SQUARED, BUT IT IS NOT VERY IMPORTANT USING GENERIC HIGH ORDER
end
%for g=1:length(svr.ParSVR.regFracOrd);EHE=EHE+svr.ParSVR.ti(g)*filtering(filtering(svr.x,svr.F{g}),conj(svr.F{g}));end
%for g=1:length(svr.ParSVR.regFracOrd);EHE=EHE+svr.ParSVR.ti(g)*filtering(svr.x,svr.F{g});end

if svr.ParSVR.tiSh~=0 && ~initStep;EHE=EHE+(svr.ParSVR.tiSh/2)*real(shearletTransformGPUMemory(shearletTransformGPUMemory(svr.x,svr.Sh),svr.Sh,-1,svr.We,2));end

r=r-EHE;
if isfield(svr,'Precond');z=bsxfun(@times,r,svr.Precond);else z=r;end
    
p=z;
zr=real(conj(z).*r);
%zr=conj(z).*r;
rsold=sum(zr(:));
if rsold<1e-10;svr.x=z;return;end  

err=inf;
for n=1:nIt   
    svr.xx=p;
    svr=svrDecode(svrEncode(svr));
    EHE=svr.xx;
    for g=1:length(svr.ParSVR.regFracOrd)
        if svr.ParSVR.ti(g)~=0;EHE=EHE+svr.ParSVR.ti(g)*real(filtering(p,svr.F{g},1));end%%%NOTE THAT THIS SHOULD BE SQUARED, BUT IT IS NOT VERY IMPORTANT USING GENERIC HIGH ORDERS
    end
    %for g=1:length(svr.ParSVR.regFracOrd);EHE=EHE+svr.ParSVR.ti(g)*filtering(filtering(p,svr.F{g}),conj(svr.F{g}));end
    %for g=1:length(svr.ParSVR.regFracOrd);EHE=EHE+svr.ParSVR.ti(g)*filtering(p,svr.F{g});end
    if svr.ParSVR.tiSh~=0 && ~initStep;EHE=EHE+(svr.ParSVR.tiSh/2)*real(shearletTransformGPUMemory(shearletTransformGPUMemory(p,svr.Sh),svr.Sh,-1,svr.We,2));end

        
    g=rsold./sum(real(conj(p(:)).*EHE(:)));%Previous implementation was using conj(rsold) here 
    %g=conj(rsold)./sum(conj(p(:)).*EHE(:));%Previous implementation was using conj(rsold) here 
    xup=g*p;
    svr.x=svr.x+xup;    
    
    if norm~=2;svrWriteData(svr,sprintf('Step%02d-FinalReconstruction',svr.ParSVR.Step),svr.x,[],0.8,[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;end
          
    r=r-g*EHE;
    if isfield(svr,'Precond');z=bsxfun(@times,r,svr.Precond);else z=r;end
        
    zr=real(conj(z).*r);
    %zr=conj(z).*r;
    rs=sum(zr(:));
    d=rs/rsold;

    p=z+d*p;
    rsold=rs;

    if tolTyp==0;err=sqrt(sum(abs(zr(:))))/(Del*sqrt(sum(abs(x(:).*conj(x(:)))))+sqrt(bnorm2));
    elseif tolTyp==1;err=sqrt(sum(abs(zr(:)))/bnorm2);
    elseif tolTyp==2;err=sqrt(sum(abs(xup(:).*conj(xup(:))))/bnorm2);
    elseif tolTyp==3;err=sqrt(max(abs(zr(:)))/bnormI);
    elseif tolTyp==4;err=sqrt(max(abs(xup(:)))/bnormI);
    end 

    if err<tol || abs(rs)<1e-6;break;end
end

if svr.ParSVR.usePrecond==2 && norm~=1
    x=svr.x;
    svr=svrCG(svr,nIt,tolType,tol,NitTest,1,initStep);
    svr.x=x./(svr.x+eps);
end

if n==nIt && NitTest==1;fprintf('CG solver terminated without reaching convergence\n');end

