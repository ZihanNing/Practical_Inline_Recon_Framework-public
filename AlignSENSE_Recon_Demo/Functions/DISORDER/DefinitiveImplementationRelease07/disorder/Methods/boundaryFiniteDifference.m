function x=boundaryFiniteDifference(x,N,or,per)

%BOUNDARYFINITEDIFFERENCE   Imposes the boundary conditions of a finite
%difference operator
%   Y=BOUNDARYFINITEDIFFERENCE(X,N,OR,{PER})
%   X is the array on which to operate
%   N are the number of dimensions of the finite difference operator
%   OR is the finite difference order (just 1 and 2 are supported so far)
%   {PER} are the periodicity conditions: 0-> non-periodic / 1-> periodic. 
%   Defaults to 1.
%   X is the array with the boundary conditions imposed
%

if ~exist('per','var') || isempty(per);per=single(ones(1,N));end

assert(ismember(or,[1 2]),'Only first and second order finite differences are supported');
Np=length(per);per(Np+1:N)=1;

M=size(x);
y=repmat(x,[ones(1,length(M)) N]);
y(:)=1;

for l=1:N
    if ~per(l)
        if or==2;y=dynInd(y,{[1 M(l)],l},[l length(M)+1],0);else y=dynInd(y,{1,l},[l length(M)+1],0);end
    end
end
x=sum(y,length(M)+1);
