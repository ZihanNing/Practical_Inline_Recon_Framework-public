function et=extractFactorsSincRigidTransform(et,vA,dimS)

%EXTRACTFACTORSSINCRIGIDTRANSFORM extracts a subset of the k-space phase 
%multiplicative factors  required to apply a rigid transform based on sinc 
%interpolation
%   ET=EXTRACTFACTORSSINCRIGIDTRANSFORM(ET,VA,DIMS) 
%   * ET is the whole set of factors to apply the transform
%   * VA are indexes to the motion states of interest
%   * DIMS is the dimension on which the motion states are arranged
%   ** ET is a subset of factors to apply the transform 
%

if ~iscell(et{1})
    et{1}=dynInd(et{1},vA,dimS);
elseif max(vA)<=size(et{1}{1},dimS)
    for n=1:3;et{1}{n}=dynInd(et{1}{n},vA,dimS);end
end
for n=1:3
    if length(et)==5
        for m=[2:3 5];et{m}{n}=dynInd(et{m}{n},vA,dimS);end
    else
        for m=2:3;et{m}{n}=dynInd(et{m}{n},vA,dimS);end
    end
end

