function [P,SH]=preconditionMSAlignedSense(S,M)

%PRECONDITIONMSALIGNEDSENSE   Builds the preconditioner
%   P=PRECONDITIONMSALIGNEDSENSE(S,{M})  
%   * S are the sensitivities
%   * {M} is a regularization for inversion
%   ** P is the preconditioner
%   ** SH are the conjugated sensitivity maps
%

if nargin<2 || isempty(M);M=inf;end
P=1./(1e-12+sum(normm(S,[],4),5)+1./M.^2);
if nargout>=2;SH=conj(S);end
