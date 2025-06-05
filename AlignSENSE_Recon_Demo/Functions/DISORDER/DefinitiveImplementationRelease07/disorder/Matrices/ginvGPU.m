function X=ginvGPU(X,r,qu)

%GINVGPU   Computes the pseudo inverse of several rectangular matrices on the GPU. 
%The code is based on the methods described in VN Katsikis, D Pappas, "Fast 
%computing of the Moore-Penrose inverse matrix," Electronic Journal of Linear 
%Algebra, 17(44):637-650.
%   X=GINVGPU(X,R,{QU})
%   * X are the input matrices
%   * R are the matrices ranks
%   * B is a vector to postmultiply the computed pseudoinverse
%   * {QU} is a flag for quick construction of the inverse function.
%   Defaults to 0
%   * X are the pseudo-inverse matrices
%  

if nargin<3;qu=0;end

N=size(X);N(end+1:3)=1;
gpu=isa(X,'gpuArray');

if N(1)<=N(2);X=matfun(@ctranspose,X);end

G=emtimes(matfun(@ctranspose,X),X);
NG=size(G);NG(end+1:3)=1;
if ~qu
    Re=ones([1 NG(2:3)],'like',G);
    G=G+diagm(Re-zrGPU(Re,r,2));
else
    if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end
    Re=epsU*eye(NG(2),'like',G);
    G=bsxfun(@plus,G,Re);
end
G=(matfun(@ctranspose,G)+G)/2;%Not much impact apparently
X=matfun(@mldivide,G,matfun(@ctranspose,X));
if N(1)<=N(2);X=matfun(@ctranspose,X);end

