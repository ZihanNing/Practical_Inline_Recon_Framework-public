
function [x] = matchHist(x,y,pctAlign,useMatlabFun,deb)

%MATCHHIST transforms the arrays so that their histograms are matched.
%   [X]=MATCHHIST(X, Y, {PCTALIGN},{USEMATLABFUN},{DEB})
%   * X array to match.
%   * Y array to be matched to (reference).
%   * {PCTALIGN} input percentile landmarks that will be used for histrogram matching. For linear rescaling, use PCTALIGN=[0 100], which corresponds to rescaling respectively the minmum and maximum values.
%   * {USEMATLABFUN} to use the Matlab built-in function (Defaults to 0).
%   * {DEB} to debug and show histograms before and after matching (Defaults to 0).
%   ** X is histogram-matched array.
%
%   Yannick Brackenier 2022-01-30

if nargin<3 || isempty(pctAlign);pctAlign=linspace(5,95,4);end
if nargin<4 || isempty(useMatlabFun);useMatlabFun=0;end 
if nargin<5 || isempty(deb);deb=0;end

%%% Cast array if not of same class
if ~strcmp(class(x), class(y))
    warning('matchHist:: x and y must have same class. y converted from %s to %s. Avoid this for computational overhead.',class(y), class(x));
    y = cast(y, class(x));
end

%%% Plot histograms before matching
if deb>0
    if ~isreal(x); xPlot = abs(x);yPlot = abs(y);else xPlot=x;yPlot=y;end
    %%% Plot histograms
    nBins=100;
    h=figure('color','w');
    subplot(1,2,1);  histogram(xPlot(:),nBins)
    subplot(1,2,2);  histogram(yPlot(:),nBins)
    sgtitle('Histograms before');set(h,'color','w','Position',get(0,'ScreenSize'));

end

%%% Extract ranges to match based on percentiles (to avoid outlier influence)
Nperc = size(pctAlign,2);
rangeMatchx = zeros([1,Nperc],'like',real(x));%Make sure range is real as we are matching magnitude values
rangeMatchy = zeros([1,Nperc],'like',real(y));

for i = 1:Nperc
    rangeMatchx(i) = prctile(abs(x(:)),pctAlign(i)) ;
    rangeMatchy(i) = prctile(abs(y(:)),pctAlign(i)) ;
end

%%% Perform matching
if ~isreal(x); xPhase = angle(x);x=abs(x); y = abs(y); end%Save phase

if useMatlabFun%Internal Matlab function
    numBins=1000;%Hardcoded
    x = imhistmatchn(rescaleND(gather(x),linspace(0,1,Nperc),gather(rangeMatchx)),...%Rescale to [0 1]since required for imhistmatchn with type single 
                     rescaleND(gather(y),linspace(0,1,Nperc),gather(rangeMatchy)),...
                     numBins);
    x = rescaleND(x,rangeMatchy,linspace(0,1,Nperc));

else%Matching intensity ranges by piece-wise linearly rescaling
    x = rescaleND(x,rangeMatchy, rangeMatchx);
end

%%% Add phase
if ~isreal(x);x = abs(x).* exp(1i*xPhase);end

%%% Plot for debugging
if deb>0
    if ~isreal(x); xPlot = abs(x);yPlot = abs(y);else xPlot=x;yPlot=y;end
    %%% Plot histograms
    nBins=100;
    limPlot = [rangeMatchy(1),rangeMatchy(end)];
    h=figure('color','w');
    subplot(1,2,1);  histogram(xPlot(:),nBins,'BinLimits', limPlot)
    subplot(1,2,2);  histogram(yPlot(:),nBins,'BinLimits', limPlot)
    sgtitle('Histograms after');set(h,'color','w','Position',get(0,'ScreenSize'));
    if deb>1
        %%% Figures plotted to see range
        plot_(xPlot,limPlot,[],yPlot,[],limPlot);
    end
end

