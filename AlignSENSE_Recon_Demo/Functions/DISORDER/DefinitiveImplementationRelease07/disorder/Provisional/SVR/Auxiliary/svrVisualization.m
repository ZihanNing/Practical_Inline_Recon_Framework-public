function svrVisualization(x,y,tit,res)

%SVRVISUALIZATION   Allows to visualize SVR information at different steps
%during the reconstruction
%   SVR=SVRVISUALIZATION(X,Y,TIT,RES)
%   * X is the data to be visualized
%   * Y is some overlying information to be visualized
%   * TIT is the title
%   * RES is the resolution
%

if nargin<2;y=[];end
if nargin<3;tit=[];end
if nargin<4 || isempty(res);res=1;end

if iscell(x)
    NN=[16 10 6]*16/res;%Standard volume size
    NSl=5;
    NV=length(x);
    if ~isempty(y)
        for v=1:NV;x{v}=cat(5,x{v},y{v});end
    end    
    NC=size(x{1},5);
    z=zeros([NN NV NC],'like',real(x{1}));
    for v=1:NV;z(:,:,:,v,:)=resampling(abs(x{v}),NN,2);end
    sub=NN(3)/NSl;    
    z=dynInd(z,round(sub/2+(0:NSl-1)*sub)+1,3);
else
    NN=[16 10 6]*16/res;%Standard volume size    
    NSl=5;
    NV=6;
    if ~isempty(y);x=cat(5,x,y);end
    x=resampling(abs(x),NN,2);
    NC=size(x,5);
    sub=NN(3)/(NSl*NV);
    z=dynInd(x,round(sub/2+(0:(NSl*NV)-1)*sub)+1,3);
    z=reshape(z,[NN(1:2) NSl NV NC]);
end
z=permute(z,[2 1 3 4 5]);
vis4DArray(z,tit);
    