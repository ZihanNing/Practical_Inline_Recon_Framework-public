

function [Anew] = conditionMat(Aold,C,rescalingMethod)

%CONDITIONMAT transforms a matrix into one with a specified conditioning number.
%   [ANEW]=CONDITIONMAT(AOLD,{C},{RESCALINGMETHOD})
%   * AOLD the original matrix
%   * {C} the desired conditioning number
%   * {RESCALINGMETHOD} the method to rescale the singular values. Defaults to 'linear'.
%   ** ANEW the new matrix with the specified conditioning number
%

if nargin<2 || isempty(C);C=1;end
if nargin<3 || isempty(rescalingMethod);rescalingMethod='linear';end

%%% Take svd 
[U, S, V] = svd(Aold,'econ');

%%% Rescale singular values
sv = diag(S); %singular values listed without zeros      
if strcmp(rescalingMethod,'linear')%Linear scaling
    svRescaled = rescaleND(sv,[sv(1)/C sv(1)]);
elseif strcmp(rescalingMethod,'logarithmic')%Logarithmic scaling
    svRescaled = exp(linspace(0,-log(C),length(sv)));
else
    error('conditionMat:: %s rescling method not implemented.',rescalingMethod);
end

%%% Project back
SNew = dynInd( zerosL(S), {1:length(svRescaled), 1:length(svRescaled)}, 1:2, diag(svRescaled));     
Anew = U*SNew*(V');

