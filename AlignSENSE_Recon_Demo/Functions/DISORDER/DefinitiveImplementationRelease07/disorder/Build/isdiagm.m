function isd=isdiagm(X)

%ISDIAGM   Checks whether all square matrices in an array are diagonal
%   X=ISDIAGM(X)
%   * X is the input array
%   ** ISD is a flag that indicates whether they are diagonal
%

N=size(X);N(end+1:3)=1;
isd=1;
for n=1:prod(N(3:end))
    if ~isdiag(X(:,:,n));isd=0;break;end
end
