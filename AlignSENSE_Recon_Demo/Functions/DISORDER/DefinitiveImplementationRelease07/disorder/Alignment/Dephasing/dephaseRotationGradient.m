
function x=dephaseRotationGradient(x,T,D)

%DEPHASEROTATIONGRADIENT   Computes the dephase rotation gradient
%   X=DEPHASEROTATIONGRADIENT(X,T,D)
%   * X is the dephaseRotation
%   * T are the motion parameters (NTx2)
%   * D are the DCT atoms
%   * X is the computed dephase gradient
%

x=bsxfun(@times,1i*x,D);
x=bsxfun(@times,x,T);

