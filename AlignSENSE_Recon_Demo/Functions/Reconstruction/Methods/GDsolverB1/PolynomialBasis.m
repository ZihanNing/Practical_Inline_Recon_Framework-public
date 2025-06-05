
function [B, coef_idx, f_idx] = PolynomialBasis( N , nMax , debug)

% N = size of the volume

if nargin < 2 || isempty(nMax); nMax = 3;end
if nargin < 3 || isempty(debug); debug = 0;end

gpu=(gpuDeviceCount>0 && ~blockGPU);
if gpu; N = gpuArray(N); end 
N  = single(N);%enfore type single

%%% create a 3D cartesian grid
[X, Y, Z] = linearTerms(N);

%%% Create basis function and store in columns
size_flatten = [numel(X), 1];
f_idx = [1:numel(X)].';
coef_idx = [];
B = []; %zeros(numel(f_idx) , num_coef, 'like', X);

b=0;
dimStore=4;
%%% ZERO ORDER
if nMax>0
    B = cat(dimStore,B,ones(N));b=b+1; coef_idx{b}='0';
end
%%% FIRST ORDER
if nMax>1
    B = cat(dimStore,B,X);b=b+1; coef_idx{b}='X';
    B = cat(dimStore,B,Y);b=b+1; coef_idx{b}='Y';
    B = cat(dimStore,B,Z);b=b+1; coef_idx{b}='Z';
end
%%% SECOND ORDER
if nMax>1
    B = cat(dimStore,B,X.^2);b=b+1; coef_idx{b}='X^2';
    B = cat(dimStore,B,Y.^2);b=b+1; coef_idx{b}='Y^2';
    B = cat(dimStore,B,Z.^2);b=b+1; coef_idx{b}='Z^2';
    B = cat(dimStore,B,X.*Y);b=b+1; coef_idx{b}='XY';
    B = cat(dimStore,B,X.*Z);b=b+1; coef_idx{b}='XZ';
    B = cat(dimStore,B,Y.*Z);b=b+1; coef_idx{b}='YZ';
end

%%% THIRD ORDER
if nMax>2
    B = cat(dimStore,B,X.^3);b=b+1; coef_idx{b}='X^3';
    B = cat(dimStore,B,Y.^3);b=b+1; coef_idx{b}='Y^3';
    B = cat(dimStore,B,Z.^3);b=b+1; coef_idx{b}='Z^3';
    B = cat(dimStore,B,X.^2.*Y);b=b+1; coef_idx{b}='X^2Y';
    B = cat(dimStore,B,X.^2.*Z);b=b+1; coef_idx{b}='X^2Z';
    B = cat(dimStore,B,X.*Y.^2);b=b+1; coef_idx{b}='XY^2';
    B = cat(dimStore,B,Y.^2.*Z);b=b+1; coef_idx{b}='Y^2Z';
    B = cat(dimStore,B,X.*Z.^2);b=b+1; coef_idx{b}='XZ^2';
    B = cat(dimStore,B,Y.*Z.^2);b=b+1; coef_idx{b}='YZ^2';
    B = cat(dimStore,B,X.*Y.*Z);b=b+1; coef_idx{b}='XYZ';
end

B = resSub(B,1:3);

if nMax>3; warning('PolynomialBasis:: 4th order polynmials not implemented. Not added to basis.');end
if any(isnan(B(:))) && ~isequal(N,[1 1 1]); error('Basis function contain NaN');end
assert(all(isreal(B(:))),'Basis must be real');

end
