function y=finiteDifference(x,N,or,di,per,sp)

%FINITEDIFFERENCE   Computes finite differences for a given array
%   Y=FINITEDIFFERENCES(X,N,OR,{PER},{SP},{DI})
%   X is the array on which to operate
%   N are the number of dimensions of the finite difference operator
%   OR is the finite difference order (just 1 and 2 are supported so far)
%   {DI} indicates whether to take forward differences (1), or backward
%   differences (0, default)
%   {PER} are the periodicity conditions: 0-> non-periodic / 1-> periodic. 
%   Defaults to 1.
%   {SP} is the spacing of the input space, defaults to [1 1 1]
%   Y is the the finite difference information
%

if ~exist('di','var') || isempty(di);di=0;end
if ~exist('per','var') || isempty(per);per=single(ones(1,N));end
if ~exist('sp','var') || isempty(sp);sp=single(ones(1,N));end

assert(ismember(or,[1 2]),'Only first and second order finite differences are supported');

[Np,Ns]=parUnaFun({per,sp},@length);
per(Np+1:N)=1;sp(Ns+1:N)=1;

M=size(x);
if di && M(end)~=N;error('The size of the last dimension of the input (%d) must match the number of dimensions of the finite difference operator (%d)',M(end),N);end
if di==0;y=repmat(x,[ones(1,length(M)) N]);else y=dynInd(x,1,length(M));end%Backwards-Forwards
y(:)=0;
for l=1:N
    if or==2 
        if ~di
            y=dynInd(y,l,length(M)+1,(2*x-dynInd(x,[M(l) 1:M(l)-1],l)-dynInd(x,[2:M(l) 1],l))/sp(l)^2);
            if ~per(l);y=dynInd(y,{[1 M(l)],l},[l length(M)+1],0);end
        else
            y=y+(2*dynInd(x,l,length(M))-dynInd(x,{[M(l) 1:M(l)-1],l},[l length(M)])-dynInd(x,{[2:M(l) 1],l},[l length(M)]))/sp(l)^2;
        end
    else
        if ~di
            y=dynInd(y,n,length(M)+1,(x-dynInd(x,[M(l) 1:M(l)-1],l))/sp(l));
            if ~per(l);y=dynInd(y,{1,l},[l length(M)+1],0);end
        else
            y=y+(dynInd(x,l,length(M))-dynInd(x,{[2:M(l) 1],l},[l length(M)]))/sp(l);
        end
    end
end
