function X=invluGPU(X,Y)

%INVLUGPU   Performs matrix inversion using a LU decomposition
%   X=INVLUGPU(X,{Y})
%   * X is the set of input matrices
%   * {Y} are the coefficients of a system of equations that needs
%   to be solved
%   * X are the set of inverses or solutions of the system
%

if nargin<2;Y=[];end

X=luGPU(X);
%This second bit takes too much, but anyhow the first bit takes more than matlab's inv
X=triBackGPU(X,triForwGPU(X,Y));
