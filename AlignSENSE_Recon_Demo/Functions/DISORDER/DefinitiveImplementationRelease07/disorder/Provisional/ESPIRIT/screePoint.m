function [y,indy]=screePoint(x, gradTresh)

%SCREEPOINT   Finds the point of maximum curvature for spectral truncation.
%It is based on the point of ramp 1 in a virtual CDF, quite ad-hoc
%   Y=SCREEPOINT(X,{GRADTHESH})
%   * X are singular values or eigenvalues decreasingly sorted
%   * {GRADTHESH} if the threshold to use for gradient of the cumulative distribution
%   ** Y is the truncation
%   ** INDY is the index corresponding to the truncation
%

if nargin <2 || isempty(gradTresh);gradTresh=1;end

N=length(x);
F=cumsum(x);
F=F/F(end);%YB: normalise
indy=find(N*gradient(F)>gradTresh,1,'last');
y=x(indy);

