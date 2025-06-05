function [x,p]=dephaseRotation(T,D,TE,returnField)

%DEPHASEROTATION   Computes the dephasing
%   X=DEPHASEROTATION(T,C,X)
%   * T are the motion parameters (NTx2)
%   * C are the parameters (DCT domain), 2xNA with NA the number of atoms
%   * X are the DCT atoms
%   * X is the computed phase
%   * Y are the derivative fields
%

if nargin<3 || isempty(TE); TE=1;end
if nargin<4 || isempty(returnField); returnField=0;end

ndT=numDims(T);
    
x=sum(bsxfun(@times,D,convertRotation(T,'rad','deg')),ndT); % YB: D in Hz/degree

if returnField; return; end %In Hz
    
if nargout>1; p = 2*pi*TE * x;end %phase in radians
x=exp(+1i *2*pi*TE * x);% YB: 2piTE since x in Hz and now in radians

