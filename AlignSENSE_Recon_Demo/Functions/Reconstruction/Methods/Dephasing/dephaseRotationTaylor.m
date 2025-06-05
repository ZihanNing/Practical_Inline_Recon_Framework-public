
function [x,p]=dephaseRotationTaylor(T,fun,D,TE,returnB0)

%DEPHASEROTATIONTAYLOR   Computes the dephasing caused by the
%pose-dependent Taylor model.
%   X=DEPHASEROTATIONTAYLOR(T,FUN,D,{TE},{RETURNB0})
%   * T are the motion parameters (NTx2)
%   * FUN is a function handle to extract the relevant motion parameters for each Taylor map.
%   * D are the Taylor maps.
%   * {TE} is the echo time
%   * {RETURNB0} a flag to return the B0 fields in Hz instead of the phase
%   ** X is the complex array containing the B0 induced phase
%   ** P is the unwrapped phase
%

if nargin<4 || isempty(TE); TE=1;end
if nargin<5 || isempty(returnB0); returnB0=0;end

ndT=numDims(T);
rotPar = fun(T);%Rotation parameters used in the Taylor model

x=sum(bsxfun(@times,D,rotPar),ndT); %D in Hz/radian(^n) with n Taylor order used

if returnB0; return; end %In Hz
    
if nargout>1; p = 2*pi*TE * x;end % 2piTE since x in Hz and converted to radians
x=exp(+1i *2*pi*TE * x);% 2piTE since x in Hz and converted to radians

