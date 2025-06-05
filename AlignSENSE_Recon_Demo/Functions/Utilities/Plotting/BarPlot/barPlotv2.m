
function [] = barPlotv2(data, xLabel, xTickLabel, yLabel, legendNames, Title, numFig, FontSize)

%%% BARPLOT creates a barplot of the data 

if nargin < 1 || length(size(data))<2; error('Data provided not suitable.');end
numGroups = size(data,1);
numColumns = size(data,2);

if size(data,3)<2;data = cat(3, data , zeros(size(data))) ; end %Append standard devitation if not provided
if nargin < 2 || isempty(xLabel); xLabel = '';end
if nargin < 3 || isempty(xTickLabel); xTickLabel = cell(1,numGroups);xTickLabel(:) = {''};end
if nargin < 4 || isempty(yLabel); yLabel = '';end
if nargin < 5 || isempty(legendNames); legendNames = cell(1,numColumns);legendNames(:) = {''};end
if nargin < 6 || isempty(Title); Title = '';end
if nargin < 7 || isempty(numFig); numFig = 999;end
if nargin < 8 || isempty(FontSize); FontSize = 16;end
FontSizeLegend = FontSize - 2;

%%% set Colors
co=[0.8500    0.3250    0.0980;         
    0         0.4470    0.7410;    
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840
    0.4470    0.7410    0
    0.7410    0         0.4470
    0         0         0];

idxList = 0:800; % Large enough if numbers of levels increases
columnWidth = 1;% Width of an individual bar
groupSpacing = 1;% Distance between center of neighbouring groups
groupInterleave = numColumns*columnWidth+ groupSpacing;% Distance between edges of neighbouring groups

h = figure(numFig);hold on
set(h,'Name','','color', 'w');
p= {};
listPlots = zeros (1, numColumns);

for i =1:numColumns
        idxPlot =  idxList(1:numGroups)*(numColumns*columnWidth + groupSpacing) ... %First element of each group
                +  i ;  % Voor 6 segments

        mean = data (:,i,1);
        stdv = data(:,i,2); 
        
        p{i} = bar(idxPlot,mean,1/numColumns, 'FaceColor',co(i,:));
        listPlots(i) = p{i};
        
        errorbar(idxPlot, mean,stdv,'LineStyle','none','Color',[0 0 0]);
end
hold off

%%% Title
title( Title, 'interpreter','latex','FontSize', FontSize)

%%% Axes
xTicks = idxList(1:numGroups)*(numColumns*columnWidth + groupSpacing) + (numColumns*columnWidth+1)/2 ;

set(gca,'XTick',xTicks,'XTickLabel',xTickLabel);
set(gca, 'FontSize',FontSize, 'TickLabelInterpreter','latex');
ylabel(yLabel, 'interpreter', 'latex', 'FontSize', FontSize)
xlabel(xLabel, 'interpreter', 'latex', 'FontSize', FontSize)

%%% Legend
AX=legend(listPlots,legendNames,'Location','NorthEast','FontSize',FontSizeLegend,'interpreter','latex');    
%LEG = findobj(AX,'type','text');
set(AX,'Position',get(AX,'Position')+[0 0 0 0])%To the right and to the top         

%%% Set ranges
rangeLow = data(:,:,1)-data(:,:,2);
rangeUp = data(:,:,1)+data(:,:,2);
range = gather(1.1* [ min(min(rangeLow(:)),0)  max(rangeUp(:)) ]    );

xlim([0 numGroups*(numColumns*columnWidth + groupSpacing) ])
ylim([range(1) range(2)])  


                

