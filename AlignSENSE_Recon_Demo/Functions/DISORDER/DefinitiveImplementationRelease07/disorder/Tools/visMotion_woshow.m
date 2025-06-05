function visMotion_woshow(rec,voxSiz,t,folderName,fileName,yl,xl,op,numFig,titleName)

if nargin<3; t=[]; end
if nargin<4; folderName=[]; end
if nargin<5; fileName=[]; end
if nargin<6; yl=[]; end
if nargin<7; xl=[]; end
if nargin<8; op=[]; end
if nargin<9; numFig=[]; end
if nargin<10; titleName=''; end

co = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
       0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];
FontSizeA = 18; ModLab = 10; LineWidth = 2;

if ~isstruct(rec); T = rec; rec = []; rec.T = T; end
if ~isempty(t); xlab = '$t$ (s)'; else; xlab = '$t$ (au)'; end

T = rec.T;
ND = ndims(T); NT = size(T);
NS = prod(NT(1:ND-1));
T = reshape(T, [NS NT(ND)]);
if isempty(t); t = 1:NS; end
if exist('voxSiz','var') && ~isempty(voxSiz)
    T(:,1:3) = bsxfun(@times, T(:,1:3), voxSiz);
    un = 'mm';
else
    un = 'pix';
end
T(:,4:6) = 180*T(:,4:6)/pi;

if ~isempty(numFig); h = figure(numFig); clf;
else; h = figure('Visible', 'off'); end
set(h, 'Color', 'w', 'Position', [100, 100, 1200, 800]);

yyaxis left
for n = 1:3
    plot(t, T(:,n), 'Color', co(n,:), 'LineWidth', LineWidth);
    hold on
end
xlabel(xlab, 'Interpreter', 'latex', 'FontSize', FontSizeA + ModLab)
ylabel(sprintf('Translation (%s)', un), 'Interpreter', 'latex', 'FontSize', FontSizeA + ModLab)
if ~isempty(xl); xlim(xl); else; xlim([t(1), t(end)]); end
if ~isempty(yl) && size(yl,1) >= 1; ylim(yl(1,:)); end
grid on

yyaxis right
for n = 4:6
    plot(t, T(:,n), 'Color', co(n,:), 'LineWidth', LineWidth);
    hold on
end
ylabel('Rotation (deg)', 'Interpreter', 'latex', 'FontSize', FontSizeA + ModLab)
if ~isempty(yl) && size(yl,1) >= 2; ylim(yl(2,:)); end

if isfield(rec,'Par')
    if strcmp(rec.Par.Labels.FatShiftDir,'F') || strcmp(rec.Par.Labels.FatShiftDir,'H'); dire{3} = 'FH';
    elseif strcmp(rec.Par.Labels.FatShiftDir,'A') || strcmp(rec.Par.Labels.FatShiftDir,'P'); dire{3} = 'AP';
    else; dire{3} = 'LR'; end

    if strcmp(rec.Par.Labels.FoldOverDir,'HF') || strcmp(rec.Par.Labels.FoldOverDir,'FH'); dire{2} = 'FH';
    elseif strcmp(rec.Par.Labels.FoldOverDir,'PA') || strcmp(rec.Par.Labels.FoldOverDir,'AP'); dire{2} = 'AP';
    else; dire{2} = 'LR'; end

    if (strcmp(dire{3},'FH') || strcmp(dire{3},'AP')) && (strcmp(dire{2},'FH') || strcmp(dire{2},'AP')); dire{1} = 'LR'; end
    if (strcmp(dire{3},'FH') || strcmp(dire{3},'LR')) && (strcmp(dire{2},'FH') || strcmp(dire{2},'LR')); dire{1} = 'AP'; end
    if (strcmp(dire{3},'LR') || strcmp(dire{3},'AP')) && (strcmp(dire{2},'LR') || strcmp(dire{2},'AP')); dire{1} = 'FH'; end

    for m = 1:3; dire{m+3} = dire{m}; end
    aux = dire{4}; dire{4} = dire{6}; dire{6} = aux;
    for m = 1:3; dire{m} = strcat('Tra-', dire{m}); end
    for m = 4:6; dire{m} = strcat('Rot-', dire{m}); end
    legend(dire, 'Location', 'SouthEast', 'FontSize', FontSizeA - 6)
end

if ~strcmp(titleName, '')
    title(titleName, 'Interpreter', 'latex', 'FontSize', FontSizeA + ModLab);
end

if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir'); mkdir(folderName); end
    savePath = fullfile(folderName, [fileName, '.png']);
    print(h, savePath, '-dpng', '-r300');
end

close(h);
end
