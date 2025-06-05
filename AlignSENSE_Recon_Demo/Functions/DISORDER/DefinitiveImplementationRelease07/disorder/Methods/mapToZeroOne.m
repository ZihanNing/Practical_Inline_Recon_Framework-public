function x=mapToZeroOne(x,tole)

% MAPTOZEROONE maps a given soft or hard mas so that values very close to 
%0 or 1 are mapped onto these values
%   X=MAPTOZEROONE(X,{TOLE})
%   X is the array to be mapped
%   TOLE is the distance to 0-1 below which the values are mapped to 0-1
%   X is the mapped array
%

if nargin<2 || isempty(tole);tole=1e-6;end

x(x<=tole)=0;x(x>1-tole)=1;