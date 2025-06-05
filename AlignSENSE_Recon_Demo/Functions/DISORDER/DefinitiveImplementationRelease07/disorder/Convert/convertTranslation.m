function x=convertTranslation(x,sp,inp,out)

%CONVERTTRANSLATION   Converts units of translation
%   X=CONVERTTRANSLATION(X,SP,INP,OUT)
%   * X is a translation in units given by INP
%   * SP is the voxel spacing
%   * INP are the input units of translation
%   * OUT are the output units of translation
%   * X is a translation in units given by OUT
%

ND=ndims(x);
perm=1:ND;perm(2)=ND;perm(ND)=2;
sp=permute(sp,perm);
if strcmp(inp,'pix')
    if strcmp(out,'mm');x=bsxfun(@times,x,sp);end
elseif strcmp(inp,'mm')
    if strcmp(out,'pix');x=bsxfun(@times,x,1./sp);end
end