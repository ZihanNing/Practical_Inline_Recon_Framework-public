function [y,indy]=screePoint(x)

%SCREEPOINT   Finds the point of maximum curvature for spectral truncation.
%It is based on the point of ramp 1 in a virtual CDF, quite ad-hoc
%   Y=SCREEPOINT(X)
%   * X are singular values or eigenvalues decreasingly sorted
%   ** Y is the truncation
%   ** INDY is the index corresponding to the truncation
%

N=length(x);
F=cumsum(x);
F=F/F(end);
indy=find(N*gradient(F)>1,1,'last');
y=x(indy);

