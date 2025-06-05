
function [lim] = defRange(x)

%DEFAULTRANGE computes a dynamic range to use for plotting.
%   [LIM] = DEFAULTRANGE(X)
%
%   * X is the array to use to get the range from.
%   %% LIM is the dynamic range as an array of size [1x2].
%
%   Yannick Brackenier 2023-04-05

if ~isreal(x); x = abs(x);end

percLow = 5;
percUp = 95;

lim = cat(2, prctile(x(:),percLow), prctile(x(:),percUp));
lim = gather(lim);