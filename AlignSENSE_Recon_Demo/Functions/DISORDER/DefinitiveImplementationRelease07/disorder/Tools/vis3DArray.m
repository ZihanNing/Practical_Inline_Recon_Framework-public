function vis3DArray(x,pau,sub,tit)

%VIS3DARRAY   Visualizes a series of slices in a single figure potentially
%with some overlaid information
%   VIS3DARRAY(X,{PAU},{SL},{DY})
%   * X is the data for visualization
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {SUB} is the subsampling factor in the third dimension
%   * {TIT} is the title of the figure
%

if nargin<2 || isempty(pau);pau=1;end
if nargin<3;sub=[];end
if nargin<4;tit=[];end
if ~isempty(sub);x=x(:,:,floor(sub/2):sub:end,:,:);end

N=size(x);N(end+1:5)=1;
if N(5)==2;y=x(:,:,:,:,2);end%Overlay

x=x(:,:,:,:,1);
figure
set(gcf, 'Position', get(0,'Screensize'))  
if nargin>1 && ~isempty(tit) && ~verLessThan('matlab','R2018b');sgtitle(tit);else tit=[];end
c=1;
x=x-min(x(:));
x=x/max(x(:));

NCol=8;
Nres=N(1:3);
Nres(3)=floor(N(3)/NCol)*NCol;
x=resampling(x,Nres,2);
if N(5)==2;y=reshape(y,Nres,2);end
Nres(3)=floor(N(3)/NCol);
x=reshape(x,[Nres NCol]);
if N(5)==2;y=reshape(y,[Nres NCol]);end
N(3)=Nres(3);
N(4)=NCol;

for s=1:N(3)
    for t=1:N(4)
        if ~isempty(tit);subtightplot(N(3),N(4),c,[0 0],[0 0.08],[0 0]);else subtightplot(N(3),N(4),c,[0 0]);end
        %set(gca, 'LooseInset', get(gca,'TightInset'))
        %pos = get(gca, 'Position');
        %pos
        %pos(1) = 0.055;
        %pos(3) = 0.9;
        %set(gca, 'Position', pos)
        I=repmat(abs(x(:,:,s,t)),[1 1 3]);        
        if N(5)==2;I(:,:,3)=0.8*I(:,:,3)+0.2*y(:,:,s,t);end%Overlay
        imshow(I)
        c=c+1;
    end
end
if pau;pause;end