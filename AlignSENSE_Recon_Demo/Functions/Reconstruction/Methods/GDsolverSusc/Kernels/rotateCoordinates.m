
function [kx, ky, kz] = rotateCoordinates(kx, ky, kz, R, Nbackup)

if nargin < 5 || isempty(Nbackup); Nbackup=[];end

%%% SIZE OF ARRAY
N = size(kx);
if length(N)<3; N = Nbackup;end
assert(~( length(N)<3) ,'rotateCoordinates:: size of the array not correct.');

%%% TAKE OUT POSSIBLE SCALING
MS = sqrt(sum(R(1:3,1:3).^2,1));%Should be [1 1 1] for a pure rotation matrix
R = R(1:3,1:3)/diag(MS);

%%% FLATTEN
K = [kx(:)';ky(:)'; kz(:)'];

%%% ROTATE
K = inv( R )  * K;%Inverse since fetching k-space coordinates

%%% RESHAPE
kx = reshape(K(1,:).',N(1:3));
ky = reshape(K(2,:).',N(1:3));
kz = reshape(K(3,:).',N(1:3));
