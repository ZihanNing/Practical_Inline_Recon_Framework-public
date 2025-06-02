
function [xo] = flipPermute(xi,fl,perm, di)

%FLIPPERMUTE flips and permutes an array given a certain permutation and/or flipping order. 
%   [XO]=FLIPPERMUTE(XI,{FL},{PERM},{DI})
%   * XI the input array.
%   * {FL} the dimension to flip.
%   * {PERM} the permutation order.
%   * {DI} a flag to apply the forward (1) or inverse operation (0).
%   ** XO the modified array.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(fl); fl = zeros(1,ndims(xi));end
if nargin < 3 || isempty(perm); perm = 1:ndims(xi);end
if nargin < 4 || isempty(di); di = 1;end

perm(length(perm)+1:ndims(xi)) = length(perm)+1:ndims(xi);
fl(end+1:length(perm)) = 0;%In case perm has additional elements

if di==0
    %%%PERMUTE
    xi = ipermute(xi, perm);
    
    %%%FLIP
    for i =1:ndims(xi) 
        if fl(i);xi = flip(xi,i); end
    end
    
else
    %%%FLIP
    for i =1:ndims(xi) 
        if fl(i);xi = flip(xi,i); end
    end
    
    %%%PERMUTE
    xi = permute(xi, perm);
end

xo = xi;
