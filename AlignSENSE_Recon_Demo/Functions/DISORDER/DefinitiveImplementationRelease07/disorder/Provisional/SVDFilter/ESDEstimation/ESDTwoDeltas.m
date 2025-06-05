function esd=ESDTwoDeltas(Beta,N,w,eig)

%ESDTWODELTAS  computes the sample distribution of eigenvalue distributions comprised of two deltas, one of them at 1
%   ESD=ESDTWODELTAS(BETA,{N},{W},{EIG})
%   * BETA is the aspect ratio assumed to be between 0 and 1
%   equal than 1
%   * {N} is the number of points where to compute the two delta distribution
%   It defaults to 256
%   * {W} is the weight of the delta not located at 1. It defaults to 0.5
%   * {EIG} is the location of the delta that is not at 1. It defaults to 8
%   * ESD is a cell array with the estimated ESD. It contains the following fields
%       - ESD.GRID, the grid on which the MP distribution is computed
%       - ESD.DENS, the empirical spectral distribution density
%

if nargin<2 || isempty(N);N=256;end
if nargin<3 || isempty(w);w=0.5;end
if nargin<4 || isempty(eig);eig=8;end
    
assert(Beta>0 && Beta<1,'The aspect ratio should be confined to (0,1) but it is %.2f',Beta);
assert(w>0 && w<1,'The weigth of the second delta should be confined to (0,1) but it is %.2f',w);

minV=min(1,eig)*(1-sqrt(Beta)).^2;
maxV=max(1,eig)*(1+sqrt(Beta)).^2;

if numel(N)==1
    esd.grid=linspace(minV,maxV,N);
else
    esd.grid=N;
end

m=SparseStieltjes2(esd.grid',Beta,w,eig);

esd.dens=imag(m)'/pi;
