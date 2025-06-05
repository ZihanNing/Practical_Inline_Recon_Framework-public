function x=sphereDelaunay(x)

%SPHEREDELAUNAY computes a Delaunay triangulation on the unit sphere. The
%code modifies the methods in
%http://people.sc.fsu.edu/~jburkardt/m_src/sphere_delaunay/sphere_delaunay.m
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 May 2010
%
%  Author:
%
%    John Burkardt
%   X=SPHEREDELAUNAY(X)
%   * X are the points on the sphere arranged as rows. If singleton it is 
%   the number of points and the points are obtained by calling the Vogel's
%   method
%   * X are the generated faces arranged as rows whose column elements
%   index the different points in the face
%

if numel(x)==1;x=vogelSphereGrid(x);end

x=convhulln(x);
%Reverse the orientation of the triangles:
x=flip(x,2);
%Rotate so smallest entry is first
[~,indMin]=min(x,[],2);
indMin=cat(2,indMin,indMin+1,indMin+2);
indMin=mod(indMin-1,3)+1;
N=size(x,1);vN=(1:N)';
indMin=bsxfun(@plus,(indMin-1)*N,vN);
x=x(indMin);
%Sort lexically
x=lexicSort(x);

end
