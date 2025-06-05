function [y,x,E,EH,A,C,p,rs,rsold,r,err,En,Enprev,nIt]=extractReadouts(y,x,E,EH,A,C,p,rs,rsold,r,err,En,Enprev,nIt,re,rd)

%EXTRACTREADOUTS   Extracts a series of readouts for quicker reconstruction
%   [Y,X,E,EH,A,C,P,RSOLD]=EXTRACTREADOUTS(Y,X,E,EH,A,C,P,RSOLD,RE,RD)
%   * Y is the data
%   * X are the reconstructions
%   * E is an encoding operator
%   * EH is a decoding operator
%   * A is a preconditioner operator
%   * C is a masking operator
%   * P is an array in the CG algorithm
%   * RS is an array in the CG algorithm
%   * RSOLD is an array in the CG algorithm
%   * R are the residuals
%   * ERR are the errors
%   * RE are the set of readouts to be extracted as a row vector
%   * RD is the readout direction
%   * X is the reduced data
%   * X are the reduced reconstructions
%   * E is the reduced encoding operator
%   * EH is the reduced decoding operator
%   * A is the reduced preconditioner operator
%   * C is the reduced masking operator
%   * P is a reduced array in the CG algorithm
%   * RS is a reduced array in the CG algorithm
%   * RSOLD is a reduced array in the CG algorithm
%   * R are the reduced residuals
%   * ERR are the reduced errors
%

if ~isempty(x);x=dynInd(x,re,rd);end
if ~isempty(y);y=dynInd(y,re,rd);end
if ~isempty(r);r=dynInd(r,re,rd);end
E.Sf=dynInd(E.Sf,re,rd);
E.Uf{rd}.NY=length(re);E.Uf{rd}.NX=length(re);
%E.UAf{rd}=dynInd(E.UAf{rd},{re,re},1:2);
if isfield(E,'Gh')
    E.Gh.Pf=dynInd(E.Gh.Pf,re,rd);E.Gh.Af=dynInd(E.Gh.Af,re,rd);
    EH.Gh.Pb=dynInd(EH.Gh.Pb,re,rd);EH.Gh.Ab=dynInd(EH.Gh.Ab,re,rd);
end
if isfield(E,'Bf') && rd==2
    E.Bf=dynInd(E.Bf,re,rd);EH.Bb=dynInd(EH.Bb,re,rd);
end
EH.Ub{rd}.NY=length(re);EH.Ub{rd}.NX=length(re);
%EH.UAb{rd}=dynInd(EH.UAb{rd},{re,re},1:2);
A.Se=dynInd(A.Se,re,rd);
C.Ma=dynInd(C.Ma,re,rd);
p=dynInd(p,re,rd);
rs=dynInd(rs,re,rd);
rsold=dynInd(rsold,re,rd);
err=dynInd(err,re,rd);
En=dynInd(En,re,rd);
Enprev=dynInd(Enprev,re,rd);
nIt=dynInd(nIt,re,rd);
