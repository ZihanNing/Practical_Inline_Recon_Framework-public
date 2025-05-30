
function [idx] = centerIdx(N)

%CENTERIDX   Returns the index of the center of the array (in the DFT sense).
%   IDX=CENTERIDX(N)
%   * N is the array size
%   ** IDX is the index of the center
%
%   Yannick Brackenier 2022-07-17

idx = ceil((N+1)/2);