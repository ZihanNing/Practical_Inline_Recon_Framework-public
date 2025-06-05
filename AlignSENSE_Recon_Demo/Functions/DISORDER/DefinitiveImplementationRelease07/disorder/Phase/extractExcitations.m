function [y,P,x,E,EH,A,C,R,EP]=extractExcitations(y,P,x,E,EH,A,C,R,EP,sl)

%EXTRACTEXCITATIONS   Extracts a series of slices for quicker reconstruction
%   [Y,P,X,E,EH,A,C]=EXTRACTEXCITATIONS(Y,P,X,E,EH,A,C,SL)
%   * Y is the data
%   * P are the phase parameters
%   * X are the reconstructions
%   * E is an encoding operator
%   * EH is a decoding operator
%   * A is a preconditioner operator
%   * C is a masking operator
%   * R is a regularization operator
%   * EP is the phase correction structure
%   * SL are the set of excitations to be extracted as a row vector
%   * Y is the reduced data
%   * P are the reduced phase parameters
%   * X are the reduced reconstructions
%   * E is the reduced encoding operator
%   * EH is the reduced decoding operator
%   * A is the reduced preconditioner operator
%   * C is the reduced masking operator
%   * C is the reduced regularization operator
%   * EP is the reduced phase correction structure
%

MB=E.Uf{3}.NX/E.Uf{3}.NY;
vEx=sl;
vSl=bsxfun(@plus,sl,E.Uf{3}.NY*((0:MB-1)'));vSl=vSl';vSl=vSl(:)';

if ~isempty(y);y=dynInd(y,vEx,3);end
if ~isempty(P);P=dynInd(P,vSl,3);end
if ~isempty(x);x=dynInd(x,vSl,3);end
E.Sf=dynInd(E.Sf,vSl,3);
if isfield(E,'Bf');E.Bf=dynInd(E.Bf,vSl,3);end
E.Gh.Pf=dynInd(E.Gh.Pf,vSl,3);
E.Uf{3}.NY=length(sl);E.Uf{3}.NX=MB*length(sl);
E.UAf{3}=dynInd(E.UAf{3},{1,vSl},1:2);
if isfield(EH,'Bb');EH.Bb=dynInd(EH.Bb,vSl,3);end
EH.Gh.Pb=dynInd(EH.Gh.Pb,vSl,3);
EH.Ub{3}.NY=length(sl);EH.Ub{3}.NX=MB*length(sl);
EH.UAb{3}=dynInd(EH.UAb{3},{vSl,1},1:2);
A.Se=dynInd(A.Se,vSl,3);
if isfield(C,'Ma');C.Ma=dynInd(C.Ma,vSl,3);end
if isfield(R,'Ti');
    if isfield(R.Ti,'la');R.Ti.la=dynInd(R.Ti.la,vSl,3);end
end

if ~isempty(EP);
    EP.flagw=dynInd(EP.flagw,vEx,3);
    EP.cP=dynInd(EP.cP,vEx,3);  
    EP.cX=dynInd(EP.cX,vEx,3);  
    EP.En=dynInd(EP.En,vEx,3);
    EP.EnPrev=dynInd(EP.EnPrev,vEx,3);
    EP.w=dynInd(EP.w,vEx,3);
end
