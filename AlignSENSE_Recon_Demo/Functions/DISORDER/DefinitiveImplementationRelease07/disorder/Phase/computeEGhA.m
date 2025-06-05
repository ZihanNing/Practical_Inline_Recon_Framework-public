function [E,EH]=computeEGhA(E,EH)

%COMPUTEEGHA   Computes the ghosting matrix
%   [E,EH]=COMPUTEEGHA(E,EH)
%   * E is an encoding operator
%   * EH is a decoding operator
%   * E is an encoding operator
%   * EH is a decoding operator
%

NPF=size(E.Gh.Pf);
NPF(E.pe)=length(E.Gh.Mf);
E.Gh.Af=ones(NPF,'like',E.Gh.Pf);
for m=1:length(E.Gh.nD);nn=E.Gh.nD(m);
    indPEAf=(E.Gh.Mf==nn);
    E.Gh.Af=dynInd(E.Gh.Af,indPEAf,E.pe,bsxfun(@times,dynInd(E.Gh.Af,indPEAf,E.pe),dynInd(E.Gh.Pf,m,E.pe)));
end
EH.Gh.Ab=conj(E.Gh.Af);
