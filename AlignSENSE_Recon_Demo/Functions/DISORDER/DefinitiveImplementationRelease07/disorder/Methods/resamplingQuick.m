function xo=resamplingQuick(x,Nres,typ)

% RESAMPLINGQUICK resamples a given array using linear interpolation
%   X=RESAMPLINGQUICK(X,NRES)
%   * X is the array to be resampled
%   * X is the resampled array
%

if nargin<3;typ='linear';end

N=size(x);N(end+1:3)=1;
Nres(end+1:3)=1;

gpu=isa(x,'gpuArray');

%assert(length(N)==3,'Not supported image dimensionality (%d)',length(N));
assert(length(Nres)==3,'Not supported image dimensionality (%d)',length(Nres));

v=cell(1,3);
for m=1:3
    v{m}=single(linspace(1,N(m),Nres(m)));   
    if gpu;v{m}=gpuArray(v{m});end        
end

[grx,gry,grz]=meshgrid(v{2},v{1},v{3});

NX=size(x);NX(end+1:4)=1;
NO=prod(NX(4:end));
x=resSub(x,4:length(NX),NO);
xo=zeros([Nres NO],'like',x);
for n=1:NO;xo(:,:,:,n)=interp3(x(:,:,:,n),grx,gry,grz,typ);end
xo=reshape(xo,[Nres NX(4:end)]);

% 
% v=cell(1,3);u=cell(1,3);
% for m=1:3
%     v{m}=single(linspace(1,N(m),Nres(m)));   
%     u{m}=single(1:N(m));
%     if gpu;v{m}=gpuArray(v{m});u{m}=gpuArray(u{m});end        
% end
% 
% for m=1:3
%     perm=1:3;perm([m 1])=[1 m];
%     x=permute(x,perm);    
%     x=myLinearInterp(u{m},x,v{m});
%     x=permute(x,perm);    
% end
% 
% end
% 
% function x= myLinearInterp(u,x,v)%Seemed like a good idea, but it takes too long
%     [~,Bin]=histc(v,u);
%     H=diff(u);
% 
%     % Extra treatment if last element is on the boundary:
%     N=size(x);
%     if Bin(length(Bin))>=N(1);Bin(length(Bin))=N(1)-1;end
% 
%     w=Bin+(v-u(Bin))./H(Bin);
% 
%     % Interpolation parameters:
%     s=w-floor(w);
%     w=floor(w);
% 
%     % Shift frames on boundary:
%     edge=(w==N(1));
%     w(edge)=w(edge)-1;
%     s(edge)=1;           % Was: Sj(d) + 1;    
% 
%     x=bsxfun(@times,x(w,:,:),(1-s(:)))+bsxfun(@times,x(w+1,:,:),s(:));
% end
% 
% 
% 
