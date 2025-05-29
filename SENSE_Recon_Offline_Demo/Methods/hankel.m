function y=hankel(x,w,di)

%HANKEL   Maps a given k-space into a Hankel matrix
%   Y=HANKEL(X,W,{DI})
%   * X is some data to be mapped
%   * {W} is a window size for mapping
%   * {DI} is the direction for mapping, defaults to 1 (forward)
%   ** Y is the mapped data
%

if nargin<3 || isempty(di);di=1;end%OPPOSITE DIRECTION NOT IMPLEMENTED YET

ND=length(w);
N=size(x);N(end+1:ND+1)=1;
assert(all(N(1:ND)>=w),'Dimensions of data%s smaller than dimensions for calibration%s',sprintf(' %d',N(1:ND)),sprintf(' %d',w));

y=zeros([prod(N(1:ND)-w+1) prod(w) N(ND+1:end)],'like',x);
NY=size(y);NY(end+1:3)=1;

iw=1:prod(w);
sw=ind2subV(w,iw);
iX=cell(1,ND);
for s=1:ND
    iX{s}=0:N(s)-w(s);
    iX{s}=bsxfun(@plus,sw(:,s),iX{s});
end

iXi=cell(1,ND);
for n=1:prod(w)
    for s=1:ND;iXi{s}=iX{s}(n,:);end
    y=dynInd(y,n,2,reshape(dynInd(x,iXi,1:ND),[NY(1) 1 NY(3:end)]));
end

