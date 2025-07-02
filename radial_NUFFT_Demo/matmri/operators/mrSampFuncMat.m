function [AtA,Aty] = mrSampFuncMat(y,Crcvr,b0,sampTimes,phs_spha,phs_conc,phs_grid)
% Create A'*A and A'*y, where A is the dense MRI sampling matrix and y is the data vector. 
%   Modern GPUs make it feasible to store the dense array for a single slice, 
%   making non-iterative solutions possible 
%
%  [AtA,Aty] = mrSampFuncMat(data,Crcvr,b0,sampTimes,phs_spha,phs_conc,phs_grid)
%
% Outputs:
%   AtA: dense matrix that represents C'*S'*S*C, where C applies reciever
%       sensitivities and S performs MRI sampling
%   Aty: vector output of C'*S'*data, where data is the k-space data
%
% Inputs:
%   data:  kspace data
%   Crcvr: receiver sensitivity profiles
%   b0,sampTimes,phs_spha,phs_conc,phs_grid: k-space and B0 inhomogeneity
%     data. See sampHighOrder documentation for more details.
%
% (c) Corey Baron 2022
%

% Create sampling operator that does not include receivers
S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid);

% Build receiver operator
R = rcvrOp(Crcvr,0);

% Create data vector output
Aty = R'*(S'*y);
sz = size(Aty);

% Create R'*S'*S*R
R = single(R.maps);
S = reshape(single(S.kbase), numel(b0), []).';
S = exp(1i*S)/sqrt(prod(sz));
S = S'*S;
AtA = 0;
for n=1:size(R,4)
    AtA = AtA + conj(reshape(R(:,:,:,n),size(S,1),1)).*S.*reshape(R(:,:,:,n),1,size(S,2)); 
end


end