function x=resamplingU(x,Nres,fo,mirror,tr)

% RESAMPLING resamples a given array using the FFT
%   X=RESAMPLING(X,NRES,{FO},{MIRROR},{TR})
%   * X is the array to be resampled
%   * NRES is the new grid size
%   * {FO} determines whether the input/output is in k-space (1), shifted
%   k-space (2) or space (0), defaults to 0, with 3 it uses shifted k-space
%   without renormalization, with 4 uses k-space without renormalization
%   * {MIRROR} determines whether to mirror the image along a given
%   dimension, defaults to 0 for all dimensions 
%   * {TR} trims the singleton dimensions on Nres 
%   ** X is the resampled array
%   THIS IS THE RECENTLY OPTIMIZED FUNCTION, SEE RESAMPLINGOLD FOR PREVIOUS
%   FUNCTION
%

if nargin<3 || isempty(fo);fo=0;end
if nargin<5 || isempty(tr);tr=1;end

gpu=isa(x,'gpuArray');
comp=~isreal(x);

N=size(x);nDimsIn=length(N);
if tr==0;indTrim=[];else indTrim=find(Nres~=1,1,'last');end
if ~isempty(indTrim) && indTrim~=1;Nres(indTrim+1:end)=[];end
nDimsOu=length(Nres);

if nargin<4 || isempty(mirror);mirror=single(zeros(1,nDimsOu));end
if isscalar(mirror);mirror=mirror*ones(1,nDimsOu);end

nDimsIn=max(nDimsOu,nDimsIn);N(end+1:nDimsIn)=1;
Nor=N(1:nDimsOu);mirror=mirror(1:nDimsOu);mirrorin=mirror;mirrorin(:)=0;

NorM=Nor+(mirror==1).*Nor;NresM=Nres+(mirror==1).*Nres;
Nmin=min(NorM,NresM);Nmax=max(NorM,NresM);

zeroF=ceil((Nmax+1)/2);
orig=zeroF-ceil((Nmin-1)/2);
fina=zeroF+floor((Nmin-1)/2);
orig(mirror==2)=1;
fina(mirror==2)=Nmin(mirror==2);

for m=1:nDimsOu
    if Nor(m)~=Nres(m)
        %MIRRORING
        NNres=[Nres(1:m-1) NresM(m) N(m+1:end)];          
        mirrorin(m)=mirror(m);
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,1);end            

        %FOURIER TRANSFORM
        if ~fo
            if mirror(m)==2;F=dctmtx(NorM(m))/sqrt(NorM(m));else F=dftmtx(NorM(m))/NorM(m);end
            F=gather(cast(F,'like',x));
            if Nor(m)>Nres(m)
                iF=false(1,NorM(m));
                iF(orig(m):fina(m))=true;
                if mirror(m)~=2;iF=ifftshift(iF);end
                F=F(iF,:);
            end
        end
        if ismember(fo,[1 4]) && mirror(m)~=2            
             v=false(1,max(Nres(m),Nor(m)));
             v(orig(m):fina(m))=true;
             v=ifftshift(v);
             xRes=zeros(NNres,'like',x);
             if Nor(m)<Nres(m);x=dynInd(xRes,v,m,x);else x=dynInd(x,v,m);end       
        elseif fo>0
            if Nor(m)<Nres(m)
                %xRes=zeros(NNres,'like',x);
                %x=dynInd(xRes,orig(m):fina(m),m,x);
                xRes=x;
                x=zeros(NNres,'like',x);
                x=dynInd(x,orig(m):fina(m),m,xRes);
            else
                x=dynInd(x,orig(m):fina(m),m);
            end
        end
     
        %INVERSE FOURIER TRANSFORM
        if ~fo
            if mirror(m)==2;FH=(dctmtx(NresM(m)))*sqrt(NresM(m));else FH=dftmtx(NresM(m));end
            FH=gather(cast(FH,'like',x));    
            if Nor(m)<Nres(m)
                iF=false(1,NresM(m));
                iF(orig(m):fina(m))=true;
                if mirror(m)~=2;iF=ifftshift(iF);end
                FH=FH(iF,:);
            end
            F=FH'*F;   
            if ~comp;F=real(F);end
            if gpu;F=gpuArray(F);end
        end
        %APPLICATION
        if ~fo;x=aplGPU(F,x,m);
        elseif ~ismember(fo,3:4);x=x*sqrt(NresM(m)/NorM(m));%New normalization in Fourier/DCT domain, note I've seen that combination of parameters fo=1 mirror=2 does not work due to different arrangement of DCT
        end

        %INVERSE MIRRORING
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,0);end
        mirrorin(m)=0;
    end
end
if ~comp;x=real(x);end
