function [] = plotPTSignalsToCompare (pGT, p1, p2, numFig, titleName)
%just a temporary function to see how motion/PT signals differ when using
%the forward and backward model

if nargin < 4 || isempty(numFig);numFig=90;end
if nargin < 5 || isempty(titleName);titleName='';end

%%% REHAPE IF p are actual motion parameters
N = size(pGT);
isMotion = length(N)>2;

if isMotion
    %Permute so NCha x NStates
    pGT=permute(pGT, [6 5 1:4]); 
    p1=permute(p1, [6 5 1:4]); 
    p2=permute(p2, [6 5 1:4]); 
    %Make rotation angles in degrees
    pGT(4:6,:) = 180/pi*pGT(4:6,:);
    p1(4:6,:) = 180/pi*p1(4:6,:);
    p2(4:6,:) = 180/pi*p2(4:6,:);
else
    pGT=abs(pGT);
    p1=abs(p1);
    p2=abs(p2);
end
NSamples= size(pGT,2);
NCha= size(pGT,1);

%%%FIGURE
h  = figure(numFig);
clf; set(h, 'color','w');

%%% first set 
subplot(3,2,1); plot(1:NSamples, (pGT))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('True parameters')
axis tight

subplot(3,2,3); plot(1:NSamples, (p1))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Prediction - forward model')
axis tight

subplot(3,2,5); plot(1:NSamples, abs((pGT)-(p1)) )
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Difference with measurement - forward model')
axis tight

%%% second set
subplot(3,2,2); plot(1:NSamples, (pGT))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('True parameters')
axis tight

subplot(3,2,4); plot(1:NSamples, (p2))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Prediction - backward model')
axis tight

subplot(3,2,6); plot(1:NSamples, abs((pGT)-(p2)) )
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Difference with measurement - backward model')
axis tight

set(h,'color','w','Position',get(0,'ScreenSize'));
if isMotion; titleName='Motion';else;  titleName='PT';end
sgtitle(sprintf('Comparing  %s parameters for forward and backward fitting',titleName));


