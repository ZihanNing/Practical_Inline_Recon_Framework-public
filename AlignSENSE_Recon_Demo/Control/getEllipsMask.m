
function [M] = getEllipsMask(x, rTh)

if nargin<2 || isempty(rTh); rTh=[1 1 1];end

if length(rTh)<3; rTh(end+1:3)=1;end

%%% Extract only 1 volume
ND = numDims(x);
x=resSub(x,4:ND);
x=dynInd(x,1,4);

%%% Create grid
N = multDimSize(x,1:3);
pGrid = generateGrid(N,[],2);%From -1 --> 1

%%% Cartesian coordinates
cart = [];
for i=1:3
    cart{i} = pGrid{i}.*onesL(abs(x));
    cart{i} = cart{i} / rTh(i);
end%[X, Y, Z] = ndgrid(pGrid{1}(:), pGrid{2}(:), pGrid{3}(:));

%%% Spherical coordinates
[~,~,r] = cart2sph(cart{1}, cart{2}, cart{3});

%%% Compute mask
M = single(r<=1);
