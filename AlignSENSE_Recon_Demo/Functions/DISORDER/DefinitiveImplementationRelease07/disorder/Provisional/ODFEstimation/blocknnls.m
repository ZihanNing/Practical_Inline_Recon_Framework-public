function x = blocknnls(A,b,p_type)
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

if nargin == 2
    p_type = 'fixed';
end

[m,n] = size(A);

%Step 0
F = [];
G = 1:n;
x = zeros(n,1);
Atb = A'*b;
y = -Atb;
p = 10;
ninf = n + 1;
switch lower(p_type)
    case 'fixed'
        alpha = 1:n; %fixed permutation
    case 'random'
        alpha = randperm(n); %random permutation
end

noready = 1;
while noready
    %Step 1
    xF = x(F);
    yG = y(G);
    if all(xF>=0) & all(yG>=0)
        x = zeros(n,1);
        y = zeros(n,1);
        x(F) = xF(1:length(F));
        break;
    else
        H1 = F(xF < 0);
        H2 = G(yG < 0);
        H = union(H1,H2);
        if length(H) < ninf
            ninf = length(H);
        else if p >= 1
                p = p - 1;
            else %p==0
                
                switch lower(p_type)
                    case 'fixed'
                        r = max(H);
                    case 'random'
                        index = zeros(1,length(H));
                        for i=1:length(H)
                            index(i) = find(alpha == H(i));
                        end
                        r = alpha(max(index));
                end

                if ismember(r, H1)
                    H1 = r;
                    H2 = [];
                else
                    H1 = [];
                    H2 = r;
                end
            end
        end
        F = union(setdiff(F,H1),H2);
        G = union(setdiff(G,H2),H1);
    end
    %Step 2
    AF = A(:,F);
    AG = A(:,G);
    %xF = AF\b;
reg=1e-6;
    xF= (AF'*AF+reg*eye(size(AF,2)))\(AF'*b);
    yG = AG'*(AF*xF-b);
    x(F) = xF;
    y(G) = yG;
end
