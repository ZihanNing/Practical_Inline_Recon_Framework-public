% function y=QuESTmse(tau,Prob)
%
% Objective function for direct search for best fit to sample eigenvalues
%
% quadratic penalty; uses QuEST function
%
% variables = locations of population eigenvalues

function y=QuESTmse(tau,Prob)

global US_A USER_DEFINED_PARAMETERS

% compute MSE distance in sample eigenvalue space
[US_A.tausort,US_A.isort]=sort(tau);
US_A.lambda=QuEST(US_A.tausort./USER_DEFINED_PARAMETERS.tauscaling,USER_DEFINED_PARAMETERS.n);
y=mean((US_A.lambda-USER_DEFINED_PARAMETERS.lambdatarget).^2).*USER_DEFINED_PARAMETERS.objscaling;
