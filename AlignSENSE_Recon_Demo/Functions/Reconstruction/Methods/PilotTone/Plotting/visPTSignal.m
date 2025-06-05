
function [h] = visPTSignal (p, pType, NY, kIndex, plotType, outlWePT, stateSample, numFig, titleName, folderName, fileName, xLabel, yLabel, TR)

if nargin < 2 || isempty(pType);pType=[];end%If 'pTime', p has dimensions [NCha NSamples]
if nargin < 3 || isempty(NY);NY=[];end
if nargin < 4 || isempty(kIndex);kIndex=[];end
if nargin < 5 || isempty(plotType);plotType=2;end % (1) magn/phase (2) real/image (3) all
if nargin < 6 || isempty(outlWePT);outlWePT=[];end
if nargin < 7 || isempty(stateSample);stateSample=[];end 
if nargin < 8 || isempty(numFig);numFig=90;end
if nargin < 9 || isempty(titleName);titleName=[];end
if nargin < 10 || isempty(folderName);folderName=[];end
if nargin < 11 || isempty(fileName);fileName=[];end
if nargin < 12 || isempty(xLabel);xLabel='Motion states [\#]';end
if nargin < 13 || isempty(yLabel);yLabel={'S_{PT}','a.u.'};end %Variable name and unit
if nargin < 14 || isempty(TR);TR=[];end %In sec

%%% SET PARAMETERS
if isempty(pType); isSlice = size(p,4)>1;else; isSlice=strcmp(pType,'pSlice');end

FontSize=10;%Baseline fontsize
ModLab=7;%Axes labels
ModTit=7;%Subtitles
ModSupTit=10;%Suptitle
ModTick=3;%Ticks of the axes
LineWidth=1;%Linewidth of the plots
invCol=1;%1 generates white background

if plotType==3; NPlot = [4 1];else; NPlot = [2 1];end

%%% CONVERT TO P MATRIX
if isSlice
    pSlice = p;p=[];
    %%% Go from slice (in image domain) to temporal profile
    for n=1:2
        %pSlice=fftshiftGPU(pSlice,n);
        pSlice=fftGPU(pSlice,n)/(NY(n));
        pSlice=fftshiftGPU(pSlice,n);
    end

    %%% Extract temporal signal using sampling patters
    if size(kIndex,2)<2 %Is assumed to iSt
        idx = kIndex.';   
    else
        idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling %2 and 1 since permuted PE2 to first dimensino
    end

    pTime = [];
    for i=1:NY(5)%Multiple repeats
        idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
        pTime=cat(1,  pTime, ...
                      dynInd( resPop(dynInd(pSlice,i,5) ,1:2,[],1),...
                              dynInd(idx,idxToExtract,2) ,1));
    end
    p = permute( pTime, [4 1 2:3]);
    pTime=[];pSlice=[];
end
NSamples= size(p,2);
NCha= size(p,1);
time = (1:NSamples);
if ~isempty(TR); time=time*TR;end

%%% PLOT
if ~isempty(numFig);h=figure(numFig);clf;else; h = figure();end
set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

if ismember(plotType,[1 3]) %Magnitude and phase
    %Magnitude
    subplot(NPlot(1),NPlot(2),1); plot(time, abs(p),'LineWidth',LineWidth,'LineStyle','-')
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
    xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    ylabel(sprintf('$|%s|$ [%s]',yLabel{1},yLabel{2}),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    axis tight;
    title('\textbf{Magnitude}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)
    %Phase
    subplot(NPlot(1),NPlot(2),2); plot(time, angle(p),'LineWidth',LineWidth,'LineStyle','-')
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
    xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    ylabel(sprintf('\\angle$%s$ [rad]',yLabel{1}),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    axis tight;
    title('\textbf{Phase}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)
end
if plotType==2; idxToPlotRealImage =[1 2];else ;idxToPlotRealImage = [3 4];end
if plotType>1 %Real and imaginary part
    %Real
    subplot(NPlot(1),NPlot(2),idxToPlotRealImage(1)); plot(time, real(p),'LineWidth',LineWidth,'LineStyle','-')
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
    xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    ylabel(sprintf('Re($%s$) [%s]',yLabel{1},yLabel{2}),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    axis tight;
    title('\textbf{Real part}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)
    %Imaginary
    subplot(NPlot(1),NPlot(2),idxToPlotRealImage(2)); plot(time, imag(p),'LineWidth',LineWidth,'LineStyle','-')
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
    xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    ylabel(sprintf('Im($%s$) [%s]',yLabel{1},yLabel{2}),'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    axis tight;
    title('\textbf{Imaginary part}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)
end

%%% ADD BOX WHERE OUTLIERS ARE DETECTED
if ~isempty(outlWePT)
    if isempty(stateSample); stateSample=1:length(outlWePT);end
    subplot(NPlot(1),NPlot(2),1);%Only plot on first one
    for i=find(outlWePT==1)
        x = find(stateSample==i);

        factBox=1;
        boxx=[x(1)-.5  x(1)-.5 x(end)+.5 x(end)+.5]; %+1 because otherwise neighbourhing boxes have a small white space in between- only for visual purposes
        if plotType==2 %Plot outliers on real subplot
            boxy = zerosL(boxx);
            boxy([1,4])= gather(min(real(p(:)))*factBox);
            boxy([2,3])= gather(max(real(p(:)))*factBox);            
        else %Plot outliers on magnitude subplot
            boxy=[0    1    1      0     ]*max(abs(p(:)))*factBox;
        end
        patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none');%'EdgeColor',[.2 0 0])% 'LineWidth',0.00001,   
    end
end

%%% SUPTITLE
if ~strcmp(titleName,''); sgtitle(titleName,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModSupTit);end
    
%%% SAVE
if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir');mkdir(folderName);end
    saveFig(strcat(folderName,filesep,fileName), [], [], [], 0, [], [], '.png');
end


