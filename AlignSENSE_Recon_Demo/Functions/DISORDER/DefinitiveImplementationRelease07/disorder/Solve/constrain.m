function x=constrain(x,C)

%CONSTRAIN   Constrains a reconstruction solution according to C
%   X=CONSTRAIN(X,C)
%   * X is the reconstruction before applying the constrains
%   * C is the structure of constrains
%   ** X is the reconstruction solution after applying the constrains
%

if isempty(C);return;end

if isfield(C,'Ma');x=bsxfun(@times,x,C.Ma);end %YB: mask x
if isfield(C,'Fi');x(C.Fi.Ma)=C.Fi.Va(C.Fi.Ma);end%YB  Replaces values at mask C.Fi.Ma with predefined values in C.Fi.Va - not sure yet where and how they are calculated
if isfield(C,'Re') && C.Re;x=real(x);end%YB: real part
if isfield(C,'Po') && C.Po;x=max(x,0);end%YB: positive part
if isfield(C,'mV');x=bsxfun(@minus,x,multDimMea(x,C.mV));end%YB: substract mean in dimensions C.mV (e.g. used in motion constrain, as 5th dimension stores the motion states)
if isfield(C,'mD');x=bsxfun(@minus,x,multDimMed(x,C.mD));end%YB: substract median in dimensions C.mD
