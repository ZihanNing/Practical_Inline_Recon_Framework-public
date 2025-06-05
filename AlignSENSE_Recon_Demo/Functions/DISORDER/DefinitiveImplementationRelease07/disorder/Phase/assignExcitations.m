function [PF,xF,EF,EHF,EPF]=assignExcitations(PF,xF,EF,EHF,EPF,P,x,E,EH,EP,sl)

%ASSIGNEXCITATIONS   Assigns particular slice reconstructions or estimates to general arrays
%   [PF,XF,EF,EHF,EPF]=ASSIGNEXCITATIONS(PF,XF,EF,EHF,EPF,P,X,E,EH,EP,SL)
%   * PF are the phase parameters
%   * XF are the reconstructions
%   * EF is an encoding operator
%   * EHF is a decoding operator
%   * EPF is a phase correction structure
%   * P are the reduced phase parameters
%   * X are the reduced reconstructions
%   * E is a reduced encoding operator
%   * EH is a reduced decoding operator
%   * EP is a reduced phase correction structure
%   * SL are the set of excitations to be assigned as a row vector
%   * PF are the phase parameters
%   * XF are the reconstructions
%   * EF is an encoding operator
%   * EHF is a decoding operator
%   * EPF is a phase correction structure
%

MB=EF.Uf{3}.NX/EF.Uf{3}.NY;
vEx=sl;
vSl=bsxfun(@plus,sl,EF.Uf{3}.NY*((0:MB-1)'));vSl=vSl';vSl=vSl(:)';

PF=dynInd(PF,vSl,3,P);
xF=dynInd(xF,vSl,3,x);
EF.Gh.Pf=dynInd(EF.Gh.Pf,vSl,3,E.Gh.Pf);
EHF.Gh.Pb=dynInd(EHF.Gh.Pb,vSl,3,EH.Gh.Pb);

EPF.flagw=dynInd(EPF.flagw,vEx,3,EP.flagw);
EPF.cP=dynInd(EPF.cP,vEx,3,EP.cP);  
EPF.cX=dynInd(EPF.cX,vEx,3,EP.cX);  
EPF.w=dynInd(EPF.w,vEx,3,EP.w);
EPF.En=dynInd(EPF.En,vEx,3,EP.En);
EPF.EnPrev=dynInd(EPF.EnPrev,vEx,3,EP.EnPrev);
