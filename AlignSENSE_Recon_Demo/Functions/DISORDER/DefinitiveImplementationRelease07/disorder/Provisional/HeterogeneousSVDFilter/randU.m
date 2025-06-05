function U=randU(N,re,gpu)

%RANDU   Generates a random unitary matrix of size N. The code is based on
%the implementation in http://www.dr-qubit.org/matlab/randU.m
%   U=RANDU(N,{RE},{GPU})
%   * N is the size of the matrix
%   * {RE} indicates whether the matrix should be orthogonal instead of
%   unitary, it defaults to 0
%   * {GPU} indicates whether to generate a GPU array
%   ** U is the generated matrix
%
% RANDU   Generates a random unitary matrix of size N
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    U = RANDU(N) generates a random N x N unitary matrix,
%    distributed uniformly according to the Haar measure.
%

if nargin<2 || isempty(re);re=0;end
if nargin<3 || isempty(gpu);gpu=gpuDeviceCount && ~blockGPU;end

N(end+1:3)=1;
Nmax=N;
Nmax(1:2)=max(N(1:2));
U=randn(Nmax,'single');
if ~re;U=U+1i*randn(Nmax,'single');U=U/sqrt(2);end
for n=1:prod(N(3:end))
    [Q,R]=qr(U(:,:,n));
    R=diag(diag(R)./abs(diag(R)));
    U(:,:,n)=Q*R;
end
if gpu;U=gpuArray(U);end
U=dynInd(U,{1:N(1),1:N(2)},1:2);