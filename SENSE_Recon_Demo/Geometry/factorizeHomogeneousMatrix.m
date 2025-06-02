function [T,R,S]=factorizeHomogeneousMatrix(MT)

%FACTORIZEHOMOGENEOUSMATRIX  Factorizes a homogeneous matrix into
%translation rotation and scaling matrices. We have MT=T*R*S;
%   [T,R,S]=FACTORIZEHOMOGENEOUSMATRIX(MT)
%   * MT is a set of 4x4 homogeneous matrix
%   ** T is a set of translations
%   ** R is a set of rotations
%   ** S is a set of scalings
%

N=size(MT);N(end+1:3)=1;

T=eye(4,'like',MT);
T=repmat(T,[1 1 N(3:end)]);
T=dynInd(T,{1:3,4},1:2,dynInd(MT,{1:3,4},1:2));
MT=dynInd(MT,{1:3,4},1:2,0);
%Polar decomposition
[L,U,V]=svdm(MT,[],0);
L=diagm(L);
R=matfun(@mtimes,U,matfun(@ctranspose,V));
S=matfun(@mtimes,matfun(@mtimes,V,L),matfun(@ctranspose,V));
