
function [h] = visCalibrationMat(A, MS, type, numFig, supTitleName, folderName, fileName)

if nargin < 2 || isempty(MS);MS=[];end
if nargin < 3 || isempty(type);type='backward';end
if nargin < 4 || isempty(numFig);numFig=[];end
if nargin < 5 || isempty(supTitleName);supTitleName='';end
if nargin<6;folderName=[];end
if nargin<7;fileName=[];end

calibToDegrees=1;%Test for rot-ran imbalance

%%%SET PARAMETERS
if strcmp(type,'forward')
    error('Not implemented yet.')
elseif strcmp(type,'backward')
    assert(size(A,1)==6,'Calibration matrix does not represent 6 rigid motion parameters.')
    NGrid=[2 3];
    if ~isempty(MS);A = dynInd(A,1:3, 1, MS.'.*dynInd(A,1:3,1));end
    if ~calibToDegrees;A = dynInd(A,4:6, 1, 180/pi*dynInd(A,4:6,1));end
end
%Titles
xlab='(Virtual) Channel [\#]';
ylab=[];
if ~isempty(MS);ylab{1} = 'Tran/$S_{PT}$ [mm/a.u.]';else;ylab{1} = 'Tran/$S_{PT}$ [vox/a.u.]';end
ylab{2} = 'Rot/$S_{PT}$ [deg/a.u.]';
titleName = [];
titleName{1}= 'Calibration row for logical $t_x$';
titleName{2}= 'Calibration row for logical $t_y$';
titleName{3}= 'Calibration row for logical $t_z$';
titleName{4}= 'Calibration row for logical $\theta_z$';
titleName{5}= 'Calibration row for logical $\theta_y$';
titleName{6}= 'Calibration row for logical $\theta_x$';
titleName = strcat( '\textbf{',titleName,'}');
%Fontsizes
FontSizeA=10;%Baseline fontsize
ModLab=5;%Axes labels
ModTit=7;%Subtitles
ModSupTit=10;%Suptitle
ModTick=3;%Ticks of the axes
LineWidth=3;
invCol=1;%1 generates white background

%%%CREATE FIGURE
if ~isempty(numFig);h=figure(numFig);clf;else; h = figure();end
set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

%%%FILL SUBPLOTS
for i=1:size(A,1)
    %Create subplot
    gap=[0.1 0.05];
    if ~strcmp(supTitleName,''); marg_h=[0.1 0.1];else;marg_h=[0.1 0.05];end
    marg_w=[0.05 0.05];
    subtightplot(NGrid(1),NGrid(2),i,gap,marg_h,marg_w) 
    %Create data to plot
    row = dynInd(A,i,1);
    plot(1:length(row),row,'LineWidth',LineWidth,'LineStyle','-')
    axis([1 length(row) -inf inf]);
    %Set labels
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModTick,'TickLabelInterpreter','latex')
    xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    ylabel(ylab{ceil(i/3)},'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    title(titleName{i},'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModTit)
end

%%% SUPTITLE
if ~strcmp(supTitleName,''); sgtitle(supTitleName,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModSupTit);end
    
%%% SAVE
if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir');mkdir(folderName);end
    %export_fig(strcat(folderName,filesep,fileName,'.png'));
    saveFig(strcat(folderName,filesep,fileName), [], [], [], 0, [], [], '.png');
end
   

% labels = {};
% for m=-n:n; labels = cat(2,labels , sprintf('$SH_{%d,%d}$',n,m) );end
% leg = legend(labels)    ;
% set(leg,'Interpreter','latex','FontSize',FontSizeA+2);
% legend show
% axis tight
% AX=legend(dire);
% LEG = findobj(AX);
% set(LEG,'Color','none','TextColor',[1 1 1]*(1-invCol),'Location','NorthWest','FontSize',FontSizeA)     

