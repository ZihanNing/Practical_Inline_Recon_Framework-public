function x=buildGeometry(N,typ,par,gpu)

% BUILDGEOMETRY builds geometric objects as image masks
%   X=BUILDGEOMETRY(N,TYP,PAR,{GPU})
%   * N are the dimensions of the image mask to be generated
%   * TYP is the geometry type
%   * PAR are the parameters of the geometry in columns, different rows
%   specify different parameters to generate a collection of geometries
%   * {GPU} determines whether to use gpu computations 
%   ** X is the geometry data
%

if nargin<4 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
assert(length(N)<=3,'Geometric objects are defined up to 3 dimensions while %d dimensions have been specified',length(N));
N(end+1:3)=1;

if strcmp(typ,'Ellipsoid')
    %BOOK MEMORY
    %x=single(zeros([N size(par,1)]));
    par=permute(par,[3 2 4 1]);
    [rGrid{1},rGrid{2},rGrid{3}]=ndgrid(1:N(1),1:N(2),1:N(3));
    if gpu
        %x=gpuArray(x);
        for n=1:3;rGrid{n}=gpuArray(rGrid{n});end
    end
else
    error('Unrecognized geometric object');
end

if strcmp(typ,'Ellipsoid')
    %The parameters are specified as c1,c2,c3,r1,r2,r3,th1,th2,th3 (center,radious,rotation, respectively in pixels and degrees)
    %We assume the transform is x'=R(x-c) with R=R_3R_2R_1 which corresponds to inverse application of transforms in motion correction conventions
    %CENTER THE GRID
    for n=1:3;rGrid{n}=bsxfun(@minus,rGrid{n},dynInd(par,n,2));end
    %ORIENT THE GRID
    orGrid=cell(1,3);
    par=dynInd(par,7:9,2,convertRotation(dynInd(par,7:9,2),'deg','rad'));
    th1=dynInd(par,7,2);th2=dynInd(par,8,2);th3=dynInd(par,9,2);
    cth1=cos(th1);cth2=cos(th2);cth3=cos(th3);sth1=sin(th1);sth2=sin(th2);sth3=sin(th3);
    orGrid{1}=bsxfun(@times,cth3.*cth2,rGrid{1})+bsxfun(@times,sth3.*cth1+cth3.*sth2.*sth1,rGrid{2})+bsxfun(@times,sth3.*sth1-cth3.*sth2.*cth1,rGrid{3});
    orGrid{2}=bsxfun(@times,-sth3.*cth2,rGrid{1})+bsxfun(@times,cth3.*cth1-sth3.*sth2.*sth1,rGrid{2})+bsxfun(@times,cth3.*sth1+sth3.*sth2.*cth1,rGrid{3});
    orGrid{3}=bsxfun(@times,sth2,rGrid{1})+bsxfun(@times,-cth2.*sth1,rGrid{2})+bsxfun(@times,cth2.*cth1,rGrid{3});
    for n=1:3;rGrid{n}=[];end
    %NORMALIZE THE GRID
    for n=1:3
        orGrid{n}=bsxfun(@times,orGrid{n},1./dynInd(par,3+n,2));
        if n==1;x=orGrid{n}.^2;else x=x+orGrid{n}.^2;end
    end
    x=single(x<1);    
end
