function vis4DArray(x,pau,sub,tit)

%VIS4DARRAY   Visualizes a 4D array potentially with some overlaid
%information
%   VIS4DARRAY(X)
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
if N(5)==2
    y=x(:,:,:,:,2);
    %[FX,FY]=gradient(y);
    %y=FX.^2+FY.^2;
    %y=single(y>0.5);
    %y=morphFourier(y,4*ones(1,2),ones(1,2),zeros(1,2));
    %y=morphFourier(y,-2*ones(1,2),ones(1,2),zeros(1,2));
    %y=single(y>0.5);
end%Overlay

x=x(:,:,:,:,1);
figure
set(gcf, 'Position', get(0,'Screensize'))  
%if nargin>1 && ~isempty(tit) && ~verLessThan('matlab','R2018b');sgtitle(tit);else tit=[];end
c=1;
x=x-min(x(:));
x=x/max(x(:));
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
