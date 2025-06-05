
function [y] = FiniteDiff(x, dim, di, n)

%FINITEDIFF   Computes the finite derivative of an array along certain dimenision
%   [Y]=FINITEDIFF(X,DI,N,DIM)
% X : Signal to compute the derivative from
% {DI} : Direction of the derivatice. Forward (1 = default) or transpose% (0) differentiation
% {n} : Differentiation order
% {dim} : Dimensions over which to differentiate

if nargin< 2 || isempty(dim); dim = 1; end
if nargin< 3 || isempty(di); di = 1; end
if nargin< 4 || isempty(n); n = 1; end

NX = size(x);

%%% Dimensions to work on
NDIM = 1:ndims(x); 
NDIM_weight = ones(size(NDIM));

NDIM_exclude = setdiff(NDIM, dim);
NDIM_weight(NDIM_exclude) = 0;

y = zeros(NX);
for dim = NDIM
    
    if NDIM_weight(dim) > 0 
        x_ext = cat(dim, x, dynInd(x, 1:n, dim)); % make circular
        %x_ext = cat(dim, x, dynInd(x, NX(dim):-1:NX(dim)-n+1, dim)); % make circular
        
        if di ==1% Forward differentiation
            dif = diff(x_ext,n,dim);
            y = y + NDIM_weight(dim) * dif;
        
        elseif di==0 % Transpose differentiation
            x_ext = circshift( x_ext, +n, dim); %shift one to the right            
            dif = diff(x_ext,n,dim);
            y = y + (-1)^n *  NDIM_weight(dim) * dif; % transpose finite difference is the backward difference of the signal shifted -1 
        end
        
        
    end
end








    
    
