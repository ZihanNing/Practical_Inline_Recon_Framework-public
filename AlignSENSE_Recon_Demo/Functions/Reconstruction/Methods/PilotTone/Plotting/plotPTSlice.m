

function [] = plotPTSlice (pSlice, numFig, titleName, NY, kIndex, plotType)

if nargin < 2 || isempty(numFig);numFig=90;end
if nargin < 3 || isempty(titleName);titleName='';end
if nargin < 6 || isempty(plotType);plotType=1;end% (1) magn/phase (2) real/image (3) all

%%% Go from slice (in image domain) to temporal profile
% Move the PT signal to the fourier domain
for n=1:2
    %pSlice=fftshiftGPU(pSlice,n);
    pSlice=fftGPU(pSlice,n)/(NY(n));
    pSlice=fftshiftGPU(pSlice,n);
end

% Extract temporal signal using sampling patters
if size(kIndex,2)<2 %Is assumed to iSt
    idx = kIndex.';   
else
    idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling %2 and 1 since permuted PE2 to first dimensino
end
%                 timeMat = zeros([prod(NY(1:2)) 1]);
%                 timeMat(idx) = 1:length(idx);
%                 timeMat = resSub(timeMat, 1, NY(1:2)); 
%                 hitMat = single(timeMat>0);
pTime = [];
for i=1:NY(5)
    idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
    pTime=cat(1,  pTime, ...
        dynInd( resPop(dynInd(pSlice,i,5) ,1:2,[],1),...
                dynInd(idx, idxToExtract,2) ,1)...
                );
end

p = permute( pTime, [4 1 2:3]);

%%% Plot
h  = figure(numFig);
clf;
set(h, 'color','w');

NSamples= size(p,2);
NCha= size(p,1);

if plotType==3; grid = [4 1];else grid = [2 1];end

if ismember(plotType,[1 3]) 
    subplot(grid(1),grid(2),1); plot(1:NSamples, abs(p))
    xlabel('Motion states')
    ylabel('PT magnitude')
    axis tight;title('Magnitude')

    subplot(grid(1),grid(2),2); plot(1:NSamples, angle(p))
    xlabel('Motion states')
    ylabel('PT phase')
    axis tight;title('Phase')
end

if plotType==2; idxToPlotRealImage =[1 2];else ;idxToPlotRealImage = [3 4];end
if plotType>1
    subplot(grid(1),grid(2),idxToPlotRealImage(1)); plot(1:NSamples, real(p))
    xlabel('Motion states')
    ylabel('PT real part')
    axis tight;title('Real')

    subplot(grid(1),grid(2),idxToPlotRealImage(2)); plot(1:NSamples, imag(p))
    xlabel('Motion states')
    ylabel('PT imaginary part')
    axis tight;title('Imag')
end

set(h,'color','w','Position',get(0,'ScreenSize'));
sgtitle(titleName);
