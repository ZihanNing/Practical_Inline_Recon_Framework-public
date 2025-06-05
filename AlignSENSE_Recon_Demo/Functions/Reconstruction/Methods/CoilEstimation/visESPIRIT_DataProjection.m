

function [] = visESPIRIT_DataProjection(y, S, virCoilIncluded, rangeProj, pau, folderName, fileName, numFig)

%VISESPIRIT_DATAPROJECTION projects the data onto the estimated coil sensitivity subspace and exports the figure for later inspection.
%   VISESPIRIT_DATAPROJECTION(Y,S,{VIRCOILINCLUDED},{RANGEPROJ},{PAU},{FOLDERNAME},{FILENAME},{NUMFIG})
%  
%   * Y are the acquired coil images in k-space
%   * S are the estimated coil sensitivities in image space
%   * {VIRCOILINCLUDED} to indicate whether a body coil is included in the 4th dimension.
%   * {RANGEPROJ} is the range of errors to plot
%   * {PAU} serves to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%   * {NUMFIG} is the number of the figure to plot 
%

if nargin < 3 || isempty(virCoilIncluded); virCoilIncluded=0;end
if nargin < 4 || isempty(rangeProj); rangeProj=[];end
if nargin < 5 || isempty(pau); pau=0;end
if nargin < 6 || isempty(folderName); folderName=[];end
if nargin < 7 || isempty(fileName); fileName=[];end
if nargin < 8 || isempty(numFig); numFig=100;end

%%% CONVERT Y TO IMAGE SPACE
N = size(y);
for n=1:3
    y=ifftshift(y,n);   
    y=ifftGPU(y,n)*sqrt(N(n));%YB: normalised fft
end    

%S = dynInd(S,1,6);%For now, only take first sens maps (first eigenvector)
S = resampling(S,size(y));%Upsample sensitivities to compute projection

%%% EXTRACT PHYSICAL COILS
if virCoilIncluded
    y = dynInd(y, 2:size(y,4), 4);
    S = dynInd(S, 2:size(S,4), 4);
end

%%% CALCULATE PROJECTION
x=bsxfun(@times,sum(bsxfun(@times,conj(S),y),4),1./(normm(S,[],4)+1e-6));
P = bsxfun(@minus, y , bsxfun(@times,S, x)) ;
P = sqrt(normm(P,[],4));%RSOS

%%% PLOT PARAMETERS
FontSize = 15;

NMaps = size(P,6);
Text=cell(NMaps,1);
for i=1:NMaps;Text{i} = sprintf('Eigenvector %.0f',i);end

%%% PLOT FIRST IMAGE: the projected image
x=permute(x, [1:3,6,4:5]);
h=plotND([],abs(x), rangeProj, [],0,{[],2},[],Text,[],[],numFig);
title('Magnitude of Projection onto Eigenvectors', 'interpreter','latex', 'FontSize', FontSize);
colormap(sqrt(gray(256))); colorbar;
set(gcf, 'Position', get(0,'Screensize'))     

    %%% WRITE TO IMAGE
    if pau==1;pause;end
    if ~isempty(folderName) && ~isempty(fileName)
        if ~exist(folderName,'dir');mkdir(folderName);end
        export_fig(strcat(folderName,filesep,fileName,'_Projection.png'));
        %close(h);
    end

%%% PLOT SECOND IMAGE
P=permute(P, [1:3,6,4:5]);
h=plotND([],abs(P), rangeProj, [],0,{[],2},[],Text,[],[],numFig+1);
title('RSOS of residuals of Projection onto Eigenvectors', 'interpreter','latex', 'FontSize', FontSize);
colormap(sqrt(gray(256))); colorbar;
set(gcf, 'Position', get(0,'Screensize'))     

    %%% WRITE TO IMAGE
    if pau==1;pause;end
    if ~isempty(folderName) && ~isempty(fileName)
        if ~exist(folderName,'dir');mkdir(folderName);end
        export_fig(strcat(folderName,filesep,fileName,'_RSOS.png'));
        %close(h);
    end

