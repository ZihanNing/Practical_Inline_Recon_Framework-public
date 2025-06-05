% function dy=QuESTdmse(tau,Prob)
%
% Gradient function for direct search for best fit to sample eigenvalues
%
% quadratic penalty; uses QuEST and QuESTgrad functions
%
% variables = locations of population eigenvalues

function dy=QuESTdmse(tau,Prob)

global US_A USER_DEFINED_PARAMETERS

% compute gradient of mean squared distance between eigenvalues
p=length(US_A.tausort);
dysort=(2/p) ...
   .*QuESTgrad(US_A.tausort./USER_DEFINED_PARAMETERS.tauscaling,USER_DEFINED_PARAMETERS.n)' ...
   *(US_A.lambda-USER_DEFINED_PARAMETERS.lambdatarget).*USER_DEFINED_PARAMETERS.objscaling ...
   ./USER_DEFINED_PARAMETERS.tauscaling;

% reassign population eigenvalues back to their original positions
dy=NaN+zeros(p,1);
dy(US_A.isort)=dysort;