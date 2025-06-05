function v=ddomGPU(X)

%DDOMGPU   Computes the diagonal dominances of the set of matrices X and
%returns the result in L. They are computed as the maximum along the rows
%and columns of the minimum of the difference between the absolute value of
%the diagonal entries and the sum along the rows and columns of the
%absolute value of all the other elements normalized by the absolute value
%of the diagonal entry. Thus, diagonal dominances range between -inf and 1,
%with stable recursive inversion above 0.
%   V=DDOMGPU(X)
%   * X is the input matrix
%   * V are the returned diagonal dominances
%

gpu=isa(X,'gpuArray');
if gpu;epsU=eps(classUnderlying(X));else epsU=eps(class(X));end

ND=ndims(X);perm=1:ND;perm([2 1])=[1 2];
N=size(X);
assert(N(1)==N(2),'Diagonal dominances only defined for squared matrices while input is %d x %d',N(1),N(2));
X=abs(X);
Xd=diagm(X);
X=X-diagm(Xd);
X=bsxfun(@rdivide,X,Xd+epsU);
Xd(:)=1;
vr=min(bsxfun(@minus,Xd,sum(X,1)),[],2);
vc=min(bsxfun(@minus,permute(Xd,perm),sum(X,2)),[],1);
v=max(vr,vc);
