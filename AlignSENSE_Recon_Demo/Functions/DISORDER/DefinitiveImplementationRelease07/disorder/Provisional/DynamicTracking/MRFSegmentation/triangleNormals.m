function [v,c]=triangleNormals(x,y)

%TRIANGLENORMALS computes the normals of the triangles connecting a set of
%vertices. The length of the normals gives the area of the triangle
%   [V,C]=TRIANGLENORMALS(N)
%   * X are a set of vertices
%   * Y is a triangulation connecting the vertices
%   * V are the normals of the triangles
%   * C are the centroids of the triangles
%

e=cat(3,x(y(:,1),:),x(y(:,2),:),x(y(:,3),:));
ev=e-dynInd(e,[3 1 2],3);%Edge vectors
en=sqrt(sum(ev.^2,2));%Length of triangle sides
p=sum(en,3)/2;%Semiperimeter
a=sqrt(p.*prod(bsxfun(@minus,p,en),3));%Area of the triangles
v=cross(ev(:,:,1),ev(:,:,2),2);%Normals of the triangles
vn=sqrt(sum(v.^2,2));
v=bsxfun(@times,v,a./vn);%Normalized normals
if nargout>=2;c=mean(e,3);end%Centroids

