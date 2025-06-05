function [gm,gn,si,rv]=geometricFeatures(x,th)

%GEOMETRICFEATURES   Computes the geometric features of the intensity of a 
%volumetric grayscale image
%   [GM,GN,SI,RV]=GEOMETRICFEATURES(X,{TH})
%   X is the array on which to operate
%   * {TH} is the thresholding of the gradient magnitude to select
%   reliable surfaces. It defaults to 0.2, has to be between 0 and 1
%   GM is the image gradient modulus
%   GN is the gradient direction
%   SI is the shape index
%   RV is the radius of curvature
%

if nargin<2 || isempty(th);th=0.1;end

N=size(x);N(end+1:3)=1;
ND=max(ndims(x),3);
fl2d=N(3)==1;
if fl2d
    permfl2d=1:ND+1;permfl2d([3 1])=[1 3];
    x=permute(x,permfl2d);
end 

[grx{2},grx{1},grx{3}]=gradient(x);
x=[];
g=cat(ND+1,grx{:});%Gradient vector field of the image
gm=sqrt(sum(g.^2,ND+1));%Gradient magnitude
if nargout==1;return;end
gn=bsxfun(@rdivide,g,gm);%Gradient direction (normal)
g=[];
if nargout==2;return;end

gth=th*max(gm(:));

[grx{2},grx{1},grx{3}]=gradient(gn);
h=cat(ND+2,grx{:});%Gradient of the normal to the isosurfaces
perm=1:ND+2;perm=circshift(perm,[0 2]);

[h,gn]=parUnaFun({h,gn},@permute,perm);

[~,igns]=sort(abs(gn),1);
gu=zeros(size(gn),'like',gn);
igns2=dynInd(igns,2,1);
igns3=dynInd(igns,3,1);
igns=[];
gu=indDim(gu,igns2,1,1);
gu=indDim(gu,igns3,1,-indDim(gn,igns2,1)./indDim(gn,igns3,1));
igns1=[];igns2=[];igns3=[];
gum=sqrt(sum(gu.^2,1));
gu=bsxfun(@rdivide,gu,gum);%First orthogonal vector
gum=[];
gv=cross(gn,gu,1);%Second orthogonal vector (uxv=n)

T=cat(2,gu,gv);%Projection to the tangent plane 
gu=[];gv=[];

%whos
%size(T)
%size(h)
%pause

NT=size(T);NT(1:2)=2;NT(end+1:5)=1;
S=zeros(NT,'like',T);

for n=1:NT(5);S=dynInd(S,n,5,-emtimes(emtimes(matfun(@ctranspose,dynInd(T,n,5)),dynInd(h,n,5)),dynInd(T,n,5)));end
%S=-emtimes(emtimes(matfun(@ctranspose,T),h),T);%Shape operator
h=[];T=[];

S11=dynInd(S,[1 1],1:2);S22=dynInd(S,[2 2],1:2);S12=dynInd(S,[1 2],1:2);S21=dynInd(S,[2 1],1:2);
S=[];
K=S11.*S22-S12.*S21;%Gaussian curvature
S12=[];S21=[];
H=0.5*(S11+S22);%Mean curvature
S11=[];S22=[];
k{1}=H+sqrt(abs(H.^2-K));%Principal curvature k1
k{2}=H-sqrt(abs(H.^2-K));%Principal curvature k2
H=[];K=[];
si=0.5-(1/pi)*atan2(k{1}+k{2},k{1}-k{2});%Volumetric shape index. From 0 to 1: Cup 0.0, Rut 0.25, Saddle 0.5, Ridge 0.75, Cap 1.0
rv=sqrt((k{1}.^2+k{2}.^2)/2);%Volumetric curvedness
k=[];
rv=1./rv;%Volumetric radius of curvature
 
[gn,si,rv]=parUnaFun({gn,si,rv},@ipermute,perm);
%Non admissible points
si(rv>max(N(1:3)) | isnan(rv) | gm<gth)=0;
rv(rv>max(N(1:3)) | isnan(rv) | gm<gth)=0;

if fl2d;[gm,gn,si,rv]=parUnaFun({gm,gn,si,rv},@ipermute,permfl2d);end
