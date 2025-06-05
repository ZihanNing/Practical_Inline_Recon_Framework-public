function x=resampling(x,Nres,fo,mirror,tr)

% RESAMPLING resamples a given array using the FFT
%   X=RESAMPLING(X,NRES,{FO},{MIRROR},{TR})
%   * X is the array to be resampled
%   * NRES is the new grid size
%   * {FO} determines whether the input/output is in k-space (1), shifted
%   k-space (2) or space (0), defaults to 0
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
indTrim=find(Nres~=1,1,'last');
%if ~isempty(indTrim);Nres(indTrim+1:end)=[];end%Change
if tr==0;indTrim=[];else indTrim=find(Nres~=1,1,'last');end
if ~isempty(indTrim) && indTrim~=1;Nres(indTrim+1:end)=[];end
nDimsOu=length(Nres);

if nargin<4 || isempty(mirror);mirror=single(zeros(1,nDimsOu));end

if tr==1;assert(nDimsOu<=nDimsIn,'Resampling dimensionality (%d) is larger than image dimensionality (%d)',nDimsOu,nDimsIn);end
Nor=N(1:nDimsOu);mirror=mirror(1:nDimsOu);mirrorin=mirror;mirrorin(:)=0;

NorM=Nor+(mirror==1).*Nor;NresM=Nres+(mirror==1).*Nres; %YB: so mirror==2 not taken into account
Nmin=min(NorM,NresM);Nmax=max(NorM,NresM);

zeroF=ceil((Nmax+1)/2);
orig=zeroF-ceil((Nmin-1)/2); %YB: These are the elements to extract from the bigger F matrix to resample the array
fina=zeroF+floor((Nmin-1)/2);
orig(mirror==2)=1; %YB: Boundaries do not apply when doin the mirror=2, which goes for the application of DCT
fina(mirror==2)=Nmin(mirror==2);

for m=1:nDimsOu
    if Nor(m)~=Nres(m)
        %MIRRORING
        NNres=[Nres(1:m-1) NresM(m) N(m+1:end)];   %YB: this only applies when extracting (so fo==1 || fo ==2) and thus dimensions changed after each m accordingly       
        mirrorin(m)=mirror(m);
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,1);end

        %FOURIER TRANSFORM
        if ~fo
            if mirror(m)==2;F=dctmtx(NorM(m))/sqrt(NorM(m));else F=dftmtx(NorM(m))/NorM(m);end %YB: dividing to have both normalised! sqrt(norm) in dct and norm for dft
            if (gpu && ~isaUnderlying(x,'double')) || (~gpu && ~isa(x,'double'));F=single(F);end
            if gpu;F=gpuArray(F);end            
            if Nor(m)>Nres(m)
                iF=false(1,Nor(m));
                iF(orig(m):fina(m))=true;
                if mirror(m)~=2;iF=ifftshift(iF);end
                F=dynInd(F,iF,1);
                %if mirror(m)~=2;F=fftshift(F,1);end
                %F=dynInd(F,orig(m):fina(m),1);
                %if mirror(m)~=2;F=ifftshift(F,1);end
            end
        end
        %SHIFT
        if fo==1 && mirror(m)~=2;x=fftshift(x,m);end
        
        %ZERO-PADDING
        if fo>0
            xRes=zeros(NNres,'like',x);
            if Nor(m)<Nres(m);x=dynInd(xRes,orig(m):fina(m),m,x);else x=dynInd(x,orig(m):fina(m),m);end
        end
        
        %IFFTSHIFT
        if fo==1 && mirror(m)~=2;x=ifftshift(x,m);end           
        %INVERSE FOURIER TRANSFORM
        if ~fo
            if mirror(m)==2;FH=(dctmtx(NresM(m)))*sqrt(NresM(m));else FH=dftmtx(NresM(m));end
            if (gpu && ~isaUnderlying(x,'double')) || (~gpu && ~isa(x,'double'));FH=single(FH);end
            if gpu;FH=gpuArray(FH);end                             
            if Nor(m)<Nres(m)
                iF=false(1,Nres(m));
                iF(orig(m):fina(m))=true;
                if mirror(m)~=2;iF=ifftshift(iF);end
                FH=dynInd(FH,iF,1);                
                %if mirror(m)~=2;FH=fftshift(FH,1);end
                %FH=dynInd(FH,orig(m):fina(m),1);
                %if mirror(m)~=2;FH=ifftshift(FH,1);end
            end
            F=FH'*F;   
            if ~comp;F=real(F);end
        end
        %APPLICATION
        if ~fo
            S=size(x);S(end+1:max(nDimsIn+1,m+1))=1;
            if m~=1;x=reshape(x,[prod(S(1:m-1)) S(m) prod(S(m+1:nDimsIn))]);else x=x(:,:);end
            if m==1;x=F*x;
            elseif m~=nDimsIn;x=matfun(@mtimes,x,F.');
            else x=x*F.';
            end
            if m==1;S(m)=size(x,1);else S(m)=size(x,2);end            
            x=reshape(x,S);
        end

        %INVERSE MIRRORING
        if mirror(m)~=2;x=mirroring(x,mirrorin==1,0);end
        mirrorin(m)=0;
    end
end
if ~comp;x=real(x);end
