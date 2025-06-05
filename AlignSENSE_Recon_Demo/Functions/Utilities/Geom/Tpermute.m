

function [Pf , Pi] = Tpermute(perm)

%TPERMUTE converts a permutation order into the corresponding transformation matrix.
%   [T]=TPERMUTE(PERM)
%   * PERM the permutation order.
%   ** T the corresponding transformation matrix.
%
%   Yannick Brackenier 2022-01-30

idx = 1:16;
perm(length(perm)+1:4)=idx(length(perm)+1:4);

%%% Forward meaning (analogy with permute.m)
Pf = diag(ones(1,4));
Pf (:,perm) = Pf;

%%% Inverse (analogy with ipermute.m)
Pi = diag(ones(1,4));
Pi= Pi(:,perm);

assert(isequal(Pi, Pf'), 'Tpermute :: Not valid permuting transformation matrix.')