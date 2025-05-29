function et=precomputeFactorsLinearRigidTransform(rGrid,T)

%PRECOMPUTEFACTORSLINEARRIGIDTRANSFORM precomputes the output grid required
%to apply a rigid transform based on linear interpolation
%   [ET,ETG]=PRECOMPUTEFACTORSLINEARRIGIDTRANSFORM(RGRID,T,{DI},{CG},{GR}) 
%   * RGRID is a grid of points in the spatial domain
%   * T are the parameters of the transform
%   * ET are the parameters to apply the transform
%   * ETG are the factors to apply the derivative of the transform 
%

ndT=ndims(T);
et{1}=dynInd(T,1:3,ndT);
per=[1 3 2;
     2 1 3];
for m=1:3       
    th=dynInd(T,m+3,ndT);cth=cos(th);sth=sin(th);
    et{2}{m}=bsxfun(@plus,bsxfun(@times,cth,rGrid{per(1,m)}),-bsxfun(@times,sth,rGrid{per(2,m)}));
    et{3}{m}=bsxfun(@plus,bsxfun(@times,sth,rGrid{per(1,m)}),bsxfun(@times,cth,rGrid{per(2,m)}));
end
