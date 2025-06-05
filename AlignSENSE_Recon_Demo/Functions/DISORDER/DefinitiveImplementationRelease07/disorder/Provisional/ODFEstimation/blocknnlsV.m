function x = blocknnlsV(A,b)
%x = blocknnls(A,b,p_type);
%solves the linear least squares problem with nonnegative variables using the block principal pivoting algorithm in [1].
%
%Input:
%   A:      [MxN] matrix 
%   b:      [Mx1] vector
%   p_type: permutation type. Possible values: 
%           'random'    : random permutation. It takes long to find a solution.
%           'fixed'     : fixed permutation [1:N]. Default value. 
%
%Output
%   x:      solution
%
% [1] Portugal, Judice and Vicente, A comparison of block pivoting and
% interior point algorithms for linear least squares problems with
% nonnegative variables, Mathematics of Computation, 63(1994), pp. 625-643
%
%  
% 
%Uriel Roque
%2.5.2006

[M,N,O]=size(A);

gpu=isa(A,'gpuArray');

%Step 0
x=single(zeros(N,1,O));
if gpu;x=gpuArray(x);end
y=-matfun(@mtimes,pagefun(@transpose,A),b);
p=10;
reg=1e-6;
for o=1:O
    F=[];
    G=single(1:N);
    yo=y(:,:,o);
    xo=x(:,:,o);
    Ao=A(:,:,o);
    bo=b(:,:,o);
    ninf=N+1;
    noready=1;
    while noready
        %Step 1
        xF=xo(F);yG=yo(G);
        if all(xF>=0) && all(yG>=0)
            xo(:)=0;yo(:)=0;            
            xo(F) = xF(1:length(F));
            break;
        else
            H1 = F(xF < 0);
            H2 = G(yG < 0);
            H = union(H1,H2);
            if length(H) < ninf
                ninf = length(H);
            else
                if p>=1
                    p=p-1;
                else %p==0
                    r = max(H);
                    if ismember(r,H1)
                        H1=r;H2=[];
                    else
                        H1=[];H2=r;
                    end
                end
            end
            F = union(setdiff(F,H1),H2);
            G = union(setdiff(G,H2),H1);
        end
        %Step 2
        AF = Ao(:,F);
        AG = Ao(:,G);
        %xF = AF\b;
        regM=reg*single(eye(size(AF,2)));
        if gpu;regMG=gpuArray(regM);end
        xF= (AF'*AF+regMG)\(AF'*bo);
        yG = AG'*(AF*xF-bo);
        xo(F) = xF;
        yo(G) = yG;
    end
    x(:,:,o)=xo;
end
