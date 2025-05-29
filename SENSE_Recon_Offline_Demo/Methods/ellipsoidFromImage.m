function [par,x]=ellipsoidFromImage(x,voxsiz,inert)

%ELLIPSOIDFROMIMAGE computes the inertia ellipsoid of a 3D binary image.
%The code is adapted from
%https://uk.mathworks.com/matlabcentral/fileexchange/34104-image-ellipsoid-3d?requestedDomain=true
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-12-01,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.
%   PAR=ELLIPSOIDFROMIMAGE(X,{VOXSIZ},{INERT})
%   * X is a 3D array
%   * {VOXSIZ} are the spacings of the elements of the array
%   * {INERT} whether to compute the inertia ellipsoid. Defaults to 0, so the fit ellipsoid is computed 
%   instead
%   * PAR are the 9 ellipsoid parameters. The parameters are specified as 
%   c1,c2,c3,r1,r2,r3,th1,th2,th3 (center,radious,rotation, respectively 
%   in pixels and degrees). They correspond to the parameters used to
%   generate the ellipsoid in the BUILDGEOMETRY method
%   * X is a mask image with the ellipsoid fit
%

if nargin<2 || isempty(voxsiz);voxsiz=ones(1,3);end
if nargin<3 || isempty(inert);inert=0;end

gpu=isa(x,'gpuArray');

N=size(x);N(end+1:3)=1;

if ~inert   
    grTh=0.2;
    gm=geometricFeatures(x,grTh);
    indGrdMask=find(gm>grTh*max(gm(:)));   
    p=ind2subV(N,indGrdMask);
    par=ellipsoidFromPoints(p);        
else
    p=ind2subV(N,find(x));
    
    NP=size(p,1);%Number of points
    ci=mean(p,1);
    c=ci.*voxsiz;%Location of ellipsoid center
    p=bsxfun(@times,bsxfun(@minus,p,ci),voxsiz);%Centered points (better for numerical accuracy)
    covP=cov(p)/NP;%Compute the covariance matrix
    [U,S]=svd(covP);%To extract inertia axes
    r=2*sqrt(diag(S)*NP)';%Length of each semi axis
    [r,ir]=sort(r,'descend');%Sort axes from greater to lower
    U=U(ir,:);    
    if U(1,1)<0%Format U to ensure first axis points to positive x direction
        U=-U;
        U(:,3)=-U(:,3);%Keep matrix determinant positive
    end
    [c,r,U]=parUnaFun({c,r,U},@gather);
    a=SpinCalc('DCMtoEA123',U',1e-3,0);%Convert axes rotation matrix to Euler angles
    par=[c r a];%Concatenate parameters           
end

if nargout>1;x=buildGeometry(N,'Ellipsoid',par,gpu);end
