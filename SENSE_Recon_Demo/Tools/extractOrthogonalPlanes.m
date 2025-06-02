function [z,w]=extractOrthogonalPlanes(x,y,par)

%EXTRACTORTHOGONALPLANES   Gets the orthogonal planes from a volume for 
%segmentation visualization
%   EXTRACTORTHOGONALPLANES(X,{Y})
%   * X is the image
%   * {Y} is the segmentation
%   * {PAR} serves to slice
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {ORTHV} shows orthogonal views
%   * {SL} serves to select a given slice
%   * {DY} serves to select a given dynamic
%

x=x(:,:,:,1);

N=size(x);
x=resampling(x,max(N)*ones(1,3),3);
M=size(x);
if nargin<3 || isempty(par)
    if ~isempty(y)
        y=y(:,:,:,1);
        y=resampling(y,max(N)*ones(1,3),3);    
        par=ellipsoidFromImage(single(y>0.5));
    else
        par=mod(ceil(M/2)+1,M)+1;
    end
else
    par=par+ceil((M-N)/2);    
end
z=zeros([M(1) M(2) 3],'like',x);
if ~isempty(y);w=zeros([M(1) M(2) 3],'like',y);end
for n=1:3
    x=shiftdim(x,1);
    z(:,:,n)=dynInd(x,round(par(n)),3);
    if ~isempty(y)
        y=shiftdim(y,1);
        w(:,:,n)=dynInd(y,round(par(n)),3);
    end
end
z=z(:,:);
if ~isempty(y);w=w(:,:);else w=[];end
