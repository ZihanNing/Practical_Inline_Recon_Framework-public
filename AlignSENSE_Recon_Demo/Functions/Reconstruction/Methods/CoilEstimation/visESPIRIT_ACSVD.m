

function [] = visESPIRIT_ACSVD(V, D, threshIdx,vectorRange, pau, folderName, fileName, numFig)

%VISESPIRIT_ACSVD   plots the singular values/vectors of the fisrt step in the ESPIRiT calculation in solveESPIRIT.m
%   VISESPIRIT_ACSVD(V,D,{THRESHIDX},{VECTORRANGE},{PAU},{FOLDERNAME},{FILENAME},{NUMFIG})
%  
%   * V are the eigenvectors (filters) in k-space belonging to the eigenvalues in D
%   * D are the singular values (sorted in decreasing order)
%   * {THRESHIDX} is the index of the eigenvalue used for singular value thresholding
%   * {VECTORRANGE} is the range for the eigenvectors to plot
%   * {PAU} serves to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%   * {NUMFIG} is the number of the figure to plot 
%

if nargin < 3 || isempty(threshIdx); threshIdx=1;warning('visESPIRIT_ACSVD:: No index provided for thresholding singular values.');end
if nargin < 4 || isempty(vectorRange); vectorRange=[min(V(:)) max(V(:))];end
if nargin < 5 || isempty(pau); pau=0;end
if nargin < 6 || isempty(folderName); folderName=[];end
if nargin < 7 || isempty(fileName); fileName=[];end
if nargin < 8 || isempty(numFig); numFig=100;end

%%% RESHAPE FILTER
Vsize = size(V);
V = reshape(V,[Vsize(1)*Vsize(2),Vsize(3)]);

%%% PLOTTING PARAMETERS
LineWidth = 2;
FontSize = 15;

%%% PLOT
h=figure(numFig);
set(h,'color','w');
subplot(211)
    plot([1:length(D)],D,'LineWidth',2); hold on 
    plot([threshIdx,threshIdx],[0,D(1)],'r--','LineWidth',LineWidth)
    hLeg = legend({'Signular values','Threshold'}); set(hLeg,'Interpreter','latex','FontSize',FontSize);
    xlabel('Singular value [$\#$]','interpreter', 'latex');
    ylabel('Value [a.u.]','interpreter', 'latex');
    title('Singular values', 'interpreter', 'latex', 'FontSize',FontSize)
    axis tight

subplot(212)
    imagesc(abs(V), vectorRange ), colormap(gray(256));
    xlabel('Singular value [$\#$]','interpreter', 'latex');
    ylabel('Singular vector [a.u.]','interpreter', 'latex');
    title('Singular vectors', 'interpreter', 'latex', 'FontSize',FontSize)

set(gcf, 'Position', get(0,'Screensize'))     

%%% WRITE TO IMAGE
if pau==1;pause;end
if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir');mkdir(folderName);end
    export_fig(strcat(folderName,filesep,fileName,'.png'));
    close(h);
end
    