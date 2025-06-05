

function [] = plotPTSignals (p, numFig, titleName, outlWePT, stateSample)

if nargin < 2 || isempty(numFig);numFig=90;end
if nargin < 3 || isempty(titleName);titleName='';end
if nargin < 4 || isempty(outlWePT);outlWePT=[];end
if nargin < 5 || isempty(stateSample);stateSample=[];end

h  = figure(numFig);
clf;
set(h, 'color','w');

NSamples= size(p,2);
NCha= size(p,1);

subplot(2,1,1); plot(1:NSamples, real(p))
xlabel('Motion states')
ylabel('PT real component')
axis tight

subplot(2,1,2); plot(1:NSamples, imag(p))
xlabel('Motion states')
ylabel('PT imaginary component')
axis tight

set(h,'color','w','Position',get(0,'ScreenSize'));
sgtitle(titleName);

%%% ADD BOXED WHERE OOUTLIERS ARE DETECTED
if ~isempty(outlWePT)
    if isempty(stateSample); stateSample=1:length(outlWePT);end
    subplot(2,1,1);
    for i=find(outlWePT==1)
        x = find(stateSample==i);

        factBox=1;
        boxx=[x(1)-.5  x(1)-.5 x(end)+.5 x(end)+.5]; %+1 because otherwise neighbourhing boxes have a small white space in between- only for visual purposes
        boxy=[0    1    1      0     ]*max(abs(p(:)))*factBox;

        patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none');%'EdgeColor',[.2 0 0])% 'LineWidth',0.00001,   
    end
end