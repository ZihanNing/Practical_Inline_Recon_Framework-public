function [y,E,EH,A,C,R,x,g,c]=flipPE(y,E,EH,A,C,R,x,g,c)

%FLIPPE   Flips PE for quicker reconstruction
%   [Y,E,EH,A,C,R,X,G,C]=FLIPPE(Y,E,EH,A,C,R,{X},{G},{C})
%   * Y is the data
%   * E is an encoding operator
%   * EH is a decoding operator
%   * A is a preconditioner operator
%   * C is a masking operator
%   * X are the reconstructions
%   * G are the g-factors
%   * C are the chi2-factors
%   * Y is the data
%   * E is an encoding operator
%   * EH is a decoding operator
%   * A is a preconditioner operator
%   * C is a masking operator
%   * R is a regularization operator
%   * X are the reconstructions
%   * G are the g-factors
%   * C are the chi2-factors
%

if nargin<7;x=[];end
if nargin<8;g=[];end
if nargin<9;c=[];end

perm=1:12;perm(1:2)=[2 1];
if isfield(E,'Bf');E.Bf=permute(E.Bf,perm);EH.Bb=permute(EH.Bb,perm);end
E.Sf=permute(E.Sf,perm);
if isfield(E,'WW');E.WW=permute(E.WW,perm);end
if isfield(E,'Gh')
    E.Gh.Pf=permute(E.Gh.Pf,perm);EH.Gh.Pb=permute(EH.Gh.Pb,perm);
    %E.Gh.Af=permute(E.Gh.Af,perm);EH.Gh.Ab=permute(EH.Gh.Ab,perm);
end
NY=E.Uf{2}.NY;NX=E.Uf{2}.NX;
E.Uf{2}.NY=E.Uf{1}.NY;E.Uf{2}.NX=E.Uf{1}.NX;
EH.Ub{2}.NY=EH.Ub{1}.NY;EH.Ub{2}.NX=EH.Ub{1}.NX;
E.Uf{1}.NY=NY;E.Uf{1}.NX=NX;
EH.Ub{1}.NY=NY;EH.Ub{1}.NX=NX;
Af=E.UAf{2};
E.UAf{2}=E.UAf{1};E.UAf{1}=Af;
Ab=EH.UAb{2};
EH.UAb{2}=EH.UAb{1};EH.UAb{1}=Ab;
A.Se=permute(A.Se,perm);
if isfield(C,'Ma');C.Ma=permute(C.Ma,perm);end
if isfield(R,'Ti')
    if isfield(R.Ti,'la');R.Ti.la=permute(R.Ti.la,perm);end
end
y=permute(y,perm);
x=permute(x,perm);
g=permute(g,perm);
c=permute(c,perm);
