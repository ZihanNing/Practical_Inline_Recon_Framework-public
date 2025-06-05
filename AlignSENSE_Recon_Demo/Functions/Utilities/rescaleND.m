
function [y] = rescaleND (x,rrOu, rrIn, corr )

%RESCALEND   rescales an ND-array with a piece-wise linear function.
%   [Y]=RESCALEND(X, RROU, {RRIN},{CORR})
%   * X the array to rescale following Y = a*X+b for each segment in the piece-wise linear function.
%   * RROU output landmarks that RRIN landmarks are mapped to. size(RROU] = [1 Nlandmarks]. Defaults to [min(x(:)) max(x(:))].
%   * {RRIN} input landmarks that will be scaled to RROU. Defaults to  [min(x(:)) max(x(:))].
%   * {CORR} correction factor to affect range of rrIn, but not the bias (rrIn(1))
%   ** Y is the rescaled array.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(rrOu); rrOu= [ min( x(:) ), max( x(:) ) ];end  
if nargin < 3 || isempty(rrIn); rrIn= linspace( min( x(:) ), max( x(:) ) , size(rrOu,2) ) ;end  %Defaults to a linear rescaling for the total range
if nargin < 4 || isempty(corr); corr=1 ;end
[rrOu,rrIn] = parUnaFun({rrOu,rrIn},@gather);

%%% Initialise
y = zeros(size(x),'like',x);

%%% Loop over landmars to create piece-wise linear functions
for i = 1:size(rrIn,2)-1

    %%% Indices on which to apply rescaling
    if i == 1%First 
        idx =  x <rrIn(i+1);
        if size(rrIn,2)==2; idx = ones(size(idx),'like',idx);end%In case only 2 percentiles, include all
        
    elseif i == size(rrIn,2)-1%Last
        idx =  x >=rrIn(i);
    else
        idx =  x >=rrIn(i) & x <rrIn(i+1);
    end
    
    %%% Linear term
    a = bsxfun(@minus, dynInd(rrOu,i+1,2) , dynInd(rrOu,i,2));
    a = bsxfun( @rdivide, a, corr*(dynInd(rrIn,i+1,2) - dynInd(rrIn,i,2) )+eps);
    
    %%% Offset
    b = dynInd(rrOu,i,2) - bsxfun(@times, a,  dynInd(rrIn,i,2));

    %%% Rescale: a*x+b
    y(idx) = bsxfun(@times,a,x(idx)) + b;
    
end
