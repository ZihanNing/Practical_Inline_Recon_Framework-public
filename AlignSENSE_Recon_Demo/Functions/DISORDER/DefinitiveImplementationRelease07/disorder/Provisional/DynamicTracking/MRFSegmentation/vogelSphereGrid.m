function [x,y,z,w,v,c,f2,f1]=vogelSphereGrid(N,ra)

%VOGELSPHEREGRID uniformly distributes N points on the sphere following the
%Vogel's method
%   [X,Y,Z,W,V,C,F2,F1]=VOGELSPHEREGRID(N)
%   * N are the number of points to generate
%   * RA serves to randomize the creation. Defaults to 0
%   * X are the generated points
%   * Y is a triangulation of the generated points
%   * Z are the triangles a given point belongs to
%   * W are the points neighboring a given point
%   * V are the normals of the triangles with length given by their areas
%   * C are the centroids of the faces
%   * F2 are the triangles sharing 2 points with a given triangle
%   * F1 are the triangles sharing 1 point with a given triangle
%

if nargin<2 || isempty(ra);ra=0;end

if N~=1
    %POINTS
    ga=pi*(3-sqrt(5));%Golden angle
    if ra;ts=rand;else ts=0.5;end
    t=(ts:N)';
    z=2*t;
    t=t*ga;
    z=(1-1/N)*(1-z/(N-1));
    r=abs(sqrt(1-z.^2));
    x=[r.*cos(t) r.*sin(t) z];%Coordinates along second dimension, number of points along first dimension
else%We return the z direction (no randomization, we should probably incorporate it)
    x=[0 0 1];
end

%TRIANGULATION
if nargout>=2;y=sphereDelaunay(x);end

%TRIANGLES A POINT BELONGS TO
if nargout>=3
    z=cell(1,N);
    for n=1:N;z{n}=find(any(y==n,2));end
end

%POINTS NEIGHBORING A GIVEN POINT
if nargout>=4
    w=cell(1,N);
    for n=1:N
        w{n}=y(z{n},:);
        w{n}=unique(w{n}(:));
        w{n}(w{n}==n)=[];
        w{n}=sort(w{n});
    end
end

%NORMALS AND CENTROIDS OF THE TRIANGLES
if nargout>=5;[v,c]=triangleNormals(x,y);end

%TRIANGLES NEIGHBORING A GIVEN TRIANGLE (SHARING EITHER TWO POINTS OR ONE
%POINT)
if nargout>=7
    M=size(y,1);
    f2=cell(1,M);f1=cell(1,M);
    for m=1:M        
        pM=sum(ismember(y,y(m,:)),2);
        f2{m}=find(pM==2);
        f1{m}=find(pM==1);
    end
end        