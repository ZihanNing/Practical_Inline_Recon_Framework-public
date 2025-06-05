function y=finiteDifferenceComplex(x,N,or,per,sp)

%FINITEDIFFERENCECOMPLEX   Computes phase finite differences for a given array
%   Y=FINITEDIFFERENCECOMPLEX(X,N,OR,{PER},{SP},{DI})
%   X is the array on which to operate
%   N are the number of dimensions of the finite difference operator
%   OR is the finite difference order (just 1 and 2 are supported so far)
%   {PER} are the periodicity conditions: 0-> non-periodic / 1-> periodic. 
%   Defaults to 1.
%   {SP} is the spacing of the input space, defaults to [1 1 1]
%   Y is the the finite difference information
%

if ~exist('per','var') || isempty(per);per=single(ones(1,N));end
if ~exist('sp','var') || isempty(sp);sp=single(ones(1,N));end

assert(ismember(or,[1 2]),'Only first and second order finite differences are supported');

[Np,Ns]=parUnaFun({per,sp},@length);
per(Np+1:N)=1;sp(Ns+1:N)=1;

M=size(x);
y=repmat(x,[ones(1,length(M)) N]);
y(:)=0;
for l=1:N
    if or==2            
        y=dynInd(y,l,length(M)+1,angle((x.^2).*conj(dynInd(x,[M(l) 1:M(l)-1],l).*dynInd(x,[2:M(l) 1],l)))/sp(l)^2);
        if ~per(l);y=dynInd(y,{[1 M(l)],l},[l length(M)+1],0);end
    else
        y=dynInd(y,n,length(M)+1,angle(x.*conj(dynInd(x,[M(l) 1:M(l)-1],l)))/sp(l));
        if ~per(l);y=dynInd(y,{1,l},[l length(M)+1],0);end
    end
end
