
function [B, coef_idx, f_idx] = SH_basis( N , nn , take_cum, debug)

% N = size of the volume
% nn = list of orders to include in the Solid Harmonics basis. If scalar
% provided, this is the n_max and all orders from 

if nargin < 2 || isempty(nn); error('SH_basis: Orders provided are not compatible');end
if nargin < 3 || isempty(take_cum); take_cum = 1;end
if nargin < 4 || isempty(debug); debug = 0;end

if isequal(size(nn),[1 1]) && take_cum; nn = 0:nn;end

gpu=(gpuDeviceCount>0 && ~blockGPU);
if gpu; N = gpuArray(N); end 
N  = single(N);%enfore type single

%%% create a 3D cartesian grid
rangex = (1:N(1)) - ( floor(N(1)/2)+1); rangex =rangex / (N(1)/2);
rangey = (1:N(2)) - ( floor(N(2)/2)+1); rangey =rangey / (N(2)/2);
rangez = (1:N(3)) - ( floor(N(3)/2)+1); rangez =rangez / (N(3)/2);

[X, Y, Z] = ndgrid(rangex, rangey, rangez);
[az,elev,r] = cart2sph(X, Y, Z);

size_flatten = [numel(X), 1];

%domain = ones(size(X)); domain(r>1) = 0;
%domain = reshape( domain, size_flatten);

phi = az;
th = -elev + pi/2;
rho = r;

%%% Create basis function and store in columns
num_coef = 0; for n = nn; num_coef = num_coef + 2*n+1;end
if debug; fprintf('Number of basis functions: %d \n', num_coef);end

coef_idx = zeros(2,num_coef);
f_idx = [1:numel(X)].';
%f_idx(domain==0) = [];

B = zeros(numel(f_idx) , num_coef, 'like', X);
b = 1; % column index to store the basis functions

for n = nn
    if debug; fprintf('Storing basis order %d \n', n);end
    for m = -n:n
        Y = harmonicY(n, m, th, phi, rho,'type','real','phase',true,'norm',true);
        Y = reshape( Y, size_flatten);
        %Y(domain==0) = [];
        %Y = Y / norm(Y(:));
        
        B(:,b) = Y;%make single
        
        coef_idx(:,b) = [n;m];
        b = b +1;
    end
end

if any(isnan(B(:))) && ~isequal(N,[1 1 1]); error('B0 basis function contain NaN');end

end
