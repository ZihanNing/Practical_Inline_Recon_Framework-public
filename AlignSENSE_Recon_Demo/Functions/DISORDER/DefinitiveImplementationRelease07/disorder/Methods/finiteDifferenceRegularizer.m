function y=finiteDifferenceRegularizer(x,N,or,per,sp,sep,di)

%FINITEDIFFERENCEREGULARIZER   Computes the finite difference regularizer for a given array
%   Y=FINITEDIFFERENCEREGULARIZER(X,N,OR,{PER},{SP},{SEP},{FO})
%   X is the array on which to operate
%   N are the number of dimensions of the finite difference operator
%   OR is the finite difference order (just 1 and 2 are supported so far)
%   {PER} are the periodicity conditions: 0-> non-periodic / 1-> periodic. 
%   Defaults to 1.
%   {SP} is the spacing of the input space, defaults to [1 1 1]
%   {SEP} indicates whether the space should be separated (i.e., Laplacian)
%   or not separated (i.e., Hessian), defaults to not separated
%   {DI} indicates whether to take forward differences (1), backward
%   differences (0) or to concatenate both ([0 1], default). For using the
%   non-periodic boundary conditions, finite differences have to be taken
%   backwards
%   Y is the the finite difference information
%

if ~exist('per','var') || isempty(per);per=single(ones(1,N));end
if ~exist('sp','var') || isempty(sp);sp=single(ones(1,N));end
if ~exist('sep','var') || isempty(sep);sep=0;end
if ~exist('di','var') || isempty(di);di=[0 1];end

NDims=ndims(x);
assert(ismember(or,[1 2]),'Only first and second order finite differences are supported');

[Np,Ns]=parUnaFun({per,sp},@length);
per(Np+1:N)=1;sp(Ns+1:N)=1;

perm=cell(1,N);
for n=1:N;perm{n}=1:NDims;perm{n}(1)=n;perm{n}(n)=1;end
M=size(x);
if or~=1 || N~=3 || length(di)==1
    y=x;
    y(:)=0;
    for n=1:N 
        if M(n)>or
            xN=x;
            for d=1:length(di);xN=differ(xN,n,or,di(d));end
            y=y+xN;       
            if ~sep && or==2
                for m=n+1:N
                    xN=x;
                    for d=1:length(di);xN=differ(xN,n,1,di(d));xN=differ(xN,m,1,di(d));end
                    y=y+2*xN;
                end            
            end
        end            
    end
else
    for l=1:3
        v{l}=[M(l) 1:M(l)-1];
        w{l}=[2:M(l) 1]; 
        if ~per(l);v{l}(1)=1;w{l}(M(l))=M(l);end
    end
    sp=sp.^2;
    y=arrayfun(@fast1stOrd3D,x,dynInd(x,v{1},1),dynInd(x,w{1},1),sp(1),dynInd(x,v{2},2),dynInd(x,w{2},2),sp(2),dynInd(x,v{3},3),dynInd(x,w{3},3),sp(3));
end

function x=differ(x,l,or,fo)
    if or==1
        if fo==0%Backward difference
            if per(l);v=[size(x,l) 1:size(x,l)-1];else v=[1 1:size(x,l)-1];end
            x=x-dynInd(x,v,l);
        else%Forward difference         
            x=x-dynInd(x,[2:M(l) 1],l);          
        end
    else
        x=2*x-dynInd(x,[M(l) 1:M(l)-1],l)-dynInd(x,[2:M(l) 1],l);
        if ~fo && ~per(l);x=dynInd(x,[1 M(l)],l,0);end
    end
    x=x/(sp(l)^or);
end

end

function y=fast1stOrd3D(x,x11,x12,sp1,x21,x22,sp2,x31,x32,sp3)
    y=(2*x-x11-x12)/sp1+(2*x-x21-x22)/sp2+(2*x-x31-x32)/sp3;
end
