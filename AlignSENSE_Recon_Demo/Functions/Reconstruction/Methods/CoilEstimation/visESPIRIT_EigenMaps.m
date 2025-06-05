

function [] = visESPIRIT_EigenMaps(W, rangeEval, pau, folderName, fileName, numFig)

%VISESPIRIT_DATAPROJECTION   projects the data onto the estimated coil sensitivity subspace and exports the figure for later inspection.
%   VISESPIRIT_DATAPROJECTION(W,{RANGEEVAL},{PAU},{FOLDERNAME},{FILENAME},{NUMFIG})
%  
%   * W are the spatial eigenmaps for the ones calculated
%   * {RANGEEVAL} is the range of errors to plot
%   * {PAU} serves to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%   * {NUMFIG} is the number of the figure to plot 
%

if nargin < 2 || isempty(rangeEval); rangeEval=[];end
if nargin < 3 || isempty(pau); pau=0;end
if nargin < 4 || isempty(folderName); folderName=[];end
if nargin < 5 || isempty(fileName); fileName=[];end
if nargin < 6 || isempty(numFig); numFig=100;end

%%% PLOT PARAMETERS
FontSize = 20;

%%% PLOT
NMaps = size(W,5);
Text=cell(NMaps,1);
for i=1:NMaps;Text{i} = sprintf('Eigenmap %.0f',i);end

W = permute(W,[1:3 5 4]);
h = plotND([],abs(W), rangeEval, [],1,{[],2},[],Text,[],[],numFig);
title('Eigenmaps', 'interpreter','latex', 'FontSize', FontSize);

set(gcf, 'Position', get(0,'Screensize'))     

%%% WRITE TO IMAGE
if pau==1;pause;end
if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir');mkdir(folderName);end
    export_fig(strcat(folderName,filesep,fileName,'.png'));
    close(h);
end

