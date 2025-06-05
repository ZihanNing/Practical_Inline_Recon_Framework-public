function T=extendTransform(T)

%EXTENDTRANSFORM extends the within shot information of the transform
%adding one state by subdividing the high frequency samples
%   T=EXTENDTRANSFORM(T)
%   * T is the original transform or parameters related to each motion state
%   * T is the extended transform or parameters related to each motion state
%

NT=size(T);
T=dynInd(T,1,4,flip(dynInd(T,1,4),3));
T=dynInd(T,NT(3)+1,3,dynInd(T,NT(3),3));
T=dynInd(T,1,4,flip(dynInd(T,1,4),3));
