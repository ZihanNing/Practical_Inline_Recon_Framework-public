
function [h] = visCalibMatrix(A, type, numFig, supTitleName)

if nargin < 2 || isempty(type);type='backward';end
if nargin < 3 || isempty(numFig);numFig=[];end
if nargin < 4 || isempty(supTitleName);supTitleName='';end

calibToDegrees=1;%Test for rot-ran imbalance

%%%SET PARAMETERS
if strcmp(type,'forward')
    error('Not implemented yet.')
elseif strcmp(type,'backward')
    assert(size(A,1)==6,'Calibration matrix does not represent 6 rigid motion parameters.')
    NGrid=[2 3];
    if ~calibToDegrees;A = dynInd(A,4:6, 1, 180/pi*dynInd(A,4:6,1));end
end

xlab='(Virtual) Channel [\#]';
ylab=[];
ylab{1} = 'Translation/PT signal [vox/a.u.]';
ylab{2} = 'Rotation/PT signal [deg/a.u.]';
titleName = [];
titleName{1}= 'Calibration row for $t_x$';
titleName{2}= 'Calibration row for $t_y$';
titleName{3}= 'Calibration row for $t_z$';
titleName{4}= 'Calibration row for $\theta_z$';
titleName{5}= 'Calibration row for $\theta_y$';
titleName{6}= 'Calibration row for $\theta_x$';
titleName = strcat( '\textbf{',titleName,'}');

FontSizeA=10;
ModLab=3;
ModTit=7;
ModSupTit=10;
ModTick=2;
LineWidth=3;
invCol=1;

%%%CREATE FIGURE
if ~isempty(numFig);h=figure(numFig);clf;else; h = figure();end
set(h,'color','w','Position',get(0,'ScreenSize'));

%%%FILL SUBPLOTS
for i=1:size(A,1)
    %Create subplot
    gap=[0.1 0.05];
    if ~strcmp(supTitleName,''); marg_h=[0.1 0.1];else;marg_h=[0.1 0.05];end
    marg_w=[0.05 0.05];
    subtightplot(NGrid(1),NGrid(2), i,gap,marg_h,marg_w) 
    %Create data to plot
    row = dynInd(A,i,1);
    plot(1:length(row),row,'LineWidth',LineWidth,'LineStyle','-')
    axis([1 length(row) -inf inf]);

    %Set labels
    set(gca,'FontSize',FontSizeA+ModTick,'TickLabelInterpreter','latex')
    xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    ylabel(ylab{ceil(i/3)},'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    title(titleName{i},'Interpreter','latex','FontSize',FontSizeA+ModTit)

    
end

%%% TITLE
if ~strcmp(supTitleName,''); sgtitle(supTitleName,'Interpreter','latex','FontSize',FontSizeA+ModSupTit);end








    
%     
%     
%     labels = {};
%     for m=-n:n; labels = cat(2,labels , sprintf('$SH_{%d,%d}$',n,m) );end
%     leg = legend(labels)    ;
%     set(leg,'Interpreter','latex','FontSize',FontSizeA+2);
%     legend show
%     axis tight
