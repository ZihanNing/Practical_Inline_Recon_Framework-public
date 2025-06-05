

function [] = visSVDThresholding(U,S,V,nv,numFig,titleName,plotU,plotV)

    if nargin<4 || isempty(nv);nv = length(diag(S));end
    if nargin<5 || isempty(numFig);numFig = [];end
    if nargin<6 || isempty(titleName);titleName = '\textbf{Singular value thresholding}';end
    if nargin<7 || isempty(plotU);plotU = 0;end
    if nargin<8 || isempty(plotV);plotV = 0;end

    %%% PLOTTING PARAMETERS
    FontSize=10;%Baseline fontsize
    ModLab=10;%Axes labels
    ModLeg=7;%Axes labels
    ModTit=7;%Subtitles
    ModSupTit=10;%Suptitle
    ModTick=4;%Ticks of the axes
    LineWidth=3;
    invCol=1;%1 generates white background

    %%% PLOT SINGULAR VALUES
    if ~isempty(numFig);h=figure(numFig);clf;else; h = figure;end
    set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

    if size(S,2)>1; sv = diag(S);else; sv=S;end
    plot(1:length(sv),sv,'LineWidth',LineWidth); hold on 
    plot([nv,nv],[0,sv(1)],'r--','LineWidth',LineWidth)
    set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
    hLeg = legend({'Signular values','Threshold'}); set(hLeg,'Interpreter','latex','FontSize',FontSize+ModLeg);
    xlabel('Singular value [$\#$]','interpreter', 'latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    ylabel('Value [a.u.]','interpreter', 'latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
    title(titleName, 'interpreter', 'latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModSupTit)
    axis tight

    %%% PLOT U BASIS
    if plotU  
        if ~isempty(numFig);h=figure(numFig+1);clf;else; h = figure;end
        set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

        np = ceil(sqrt(size(U,2)));
        for i=1:size(U,2); subplot(np,np,i); plot(1:length(U(:,i)),real(U(:,i).'));end
    end
    
    %%% PLOT V BASIS
    if plotV  
        if ~isempty(numFig);h=figure(numFig+2);clf;else; h = figure;end
        set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

        np = ceil(sqrt(size(V,2)));
        for i=1:size(V,2); subplot(np,np,i); plot(1:length(V(:,i)),real(V(:,i).'));end
    end