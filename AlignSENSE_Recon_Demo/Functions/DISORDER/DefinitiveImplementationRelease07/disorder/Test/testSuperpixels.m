load('/home/lcg13/Work/ExperimentsNeuroimage18/x.mat');
N=size(x);
y=reshape(x,prod(N(1:3)),[]);
[U,S,V]=svdm(y');
V=reshape(V,N);
z=sum(abs(x).^2,4);
% 
% %Case 2D
% z=sum(z,3);
% z=z/max(z(:));
% z=reshape(z,[1 N(1) 1 N(2)]);
% z=repmat(z,[8 1 8 1]);
% z=reshape(z,8*N(1:2));
% %[L,numLabels]=superpixels(z,100);
% [L,numLabels]=superpixels(z,100,'Method','slic','Compactness',70);
% figure
% BW = boundarymask(L);
% imshow(imoverlay(z,BW,'cyan'))
% L=L(:);
% Lun=unique(L)';
% Lsi=sum(bsxfun(@eq,L,Lun),1);
% figure
% plot(Lsi)
% fprintf('Ratio of regions: %.2f\n',max(Lsi)/min(Lsi));

%Case 3D
z=z/max(z(:));
z=reshape(z,[1 N(1) 1 N(2) 1 N(3)]);
z=repmat(z,[8 1 8 1 8 1]);
z=reshape(z,8*N(1:3));
[L,numLabels]=superpixels3(z,100);%,'Method','slic','Compactness',50);
figure

imPlusBoundaries=zeros([8*N(1:2) 3 8*N(3)],'uint8');
for n=1:8*N(3)
  BW=boundarymask(L(:,:,n));
  % Create an RGB representation of this plane with boundary shown
  % in cyan.
  imPlusBoundaries(:,:,:,n) = imoverlay(z(:,:,n),BW,'cyan');
end
implay(imPlusBoundaries,5)

%BW = boundarymask(L);
%imshow(imoverlay(z,BW,'cyan'))
L=L(:);
Lun=unique(L)';
Lsi=sum(bsxfun(@eq,L,Lun),1);
figure
plot(Lsi)
fprintf('Ratio of regions: %.2f\n',max(Lsi)/min(Lsi));



%visReconstruction(sum(abs(x).^2,4),0)
%for n=1:10
%visReconstruction(V(:,:,:,n),0)
%end

