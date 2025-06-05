
function [h, hColorbar, X, B] = ...
    plotND(xRef, x, limPlot, idPlot, cmapInfo, sliceInfo, MT, labelText, M, usageM, numFig, titleText, colorbarText, FontSizeScaling, folderName, fileName)

%PLOTND is a wrapper function to facilitate plotting multiple volumes (stored in a 5D-array) in a single command with additional functionality 
%   (e.g. arrange views flexibly, colormaps, add text, mask, re-orient images, etc.).
%
%   [H, HCOLORBAR,X,B] = PLOTND({XREF}, X, {LIMPLOT}, {IDPLOT}, {CMAPINFO}, {SLICEINFO}, {MT}, {LABELTEXT}, {M}, {USAGEM}, {NUMFIG},{TITLETEXT},{COLORBARTEXT})
%   * {XREF} is a cell structure containing information about a reference volume to plot:
%       XREF{1} contains the volume, 
%       XREF{2} the dimension along which to add the reference. Defaults to 1 (=vertical).
%       XREF{3}  whether to use the colormap for the reference image. Defaults to 0.
%     If only XREF is provided, it can consist of a string referring to a NIFTI file, which will then be imported and called by plotND.
%   * X is a 5D array of volumes to plot. The 4th dimension will be plotted from top to bottom and the 5th dimension from left to right.
%   * {LIMPLOT} the range to display the arrays:
%       LIMPLOT(1:2) the range for the X array
%       LIMPLOT(3:4) can optionally contain the XREF range.
%   * {IDPLOT} cell array containing the index information about the slices to plot:
%       IDPLOT{1} the indices of the slices to plot. Defaults to the centre of the array. 
%       IDPLOT{2} a flag whether to use the id in the RAS frame (1) or the logical frame (0=dafault). 
%                 Set IDPLOT{2}=0 when IDPLOT is taken directly from NIFTI viewers (e.g.ITK-SNAP).
%   * {CMPAINFO} flag to use the jet colormap for plotting the array X (defaults to 0). Can be a:
%       Scalar (0=gray and 1 is jet colormap).
%       Char containing already existing colormaps (e.g. 'jet', 'pink', etc.).
%       Matrix of size [Ncolours 3] defining the colormap
%       Cell array containing the colors to use to form a colormap in formar HEX or NAME convention (e.g. {'red','white','blue'}).
%   * {SLICEINFO} a cell array containing information about how to plot the slices:
%       SLICEINFO{1} contains the slices to plot (sag=1,cor=2 or transverse=3). Defaults to all ([1,2,3]).
%       SLICEINFO{2} along which dimension to plot them (defaults to 2=horizontal).
%   * {MT} the orientation information in the form of the s_form from the NIFTI header (in RAS frame) to have correct orientation for the plotted arrays. If omitted, slices will not be plotted correctly.
%   * {LABELTEXT} a cell array containing strings for the titles used for the subfigures. The size of the TEXT array must be consistent with the dimensions of X.
%   * {M} the mask to use for plotting.
%   * {USAGEM} cell with USAFEM{1} a flag to indicate how to use the mask:
%       USAGEM{1}==1: Multiplicative mask.
%       USAGEM{1}==2: For plotting contours of mask with:
%             USAGEM{2} defining the contour color.
%             USAGEM{3} the contour linewidth.
%             USAGEM{4} a smoothing factor (from 0 to 1 with 1 almost no effect and 0 fully smoothed).
%       USAGEM{1}==3: Mask is used to make background black (useful when X contains RGB values).
%       USAGEM{1}==4: Mask is used to make background transparant.
%   * {NUMFIG} the number of the figure.
%   * {TITLETEXT} text to put in the title.
%   * {COLORBARTEXT} text to put in the colorbar.
%   * {FONTSIZESCALING} a 1x3 row containing the scalars to multiply with the FontSize for the (1) Title (2) Labels and (3) the colorbar.
%   ** H the handle of the figure.
%   ** HCOLORBAR the handle of the colorbar.
%   ** X the 2D array containing all the concatenated sub-figures of all the necessary slices.
%   ** B the boundary information for when the mask is used for contour plots.
%
%   Yannick Brackenier 2022-01-30

if nargin<3 || isempty(limPlot); limPlot = [];end
if nargin<4 || isempty(idPlot); idPlot = {[],0};end
if nargin<5 || isempty(cmapInfo); cmapInfo = 0;end
if nargin<6 || isempty(sliceInfo); sliceInfo = {1:3,2};end
if nargin<7 || isempty(MT); MT=[];end
if nargin<8 || isempty(labelText); labelText=[];end
if nargin<9 || isempty(M); M = [];end
if nargin<10 || isempty(usageM); usageM={0};end
if nargin<11 || isempty(numFig); numFig=[];end
if nargin<12 || isempty(titleText); titleText=[];end
if nargin<13 || isempty(colorbarText); colorbarText=[];end
if nargin<14 || isempty(FontSizeScaling); FontSizeScaling=[1 1 1];end
if nargin<15 || isempty(folderName); folderName=[];end
if nargin<16 || isempty(fileName); fileName=[];end

%% Load NIFTI file if string provided
if nargin==1 || isempty(x)%Including isempty(x) allows you use all other functionality of plotND (MT will be over-written)
    if ischar(xRef)
        if contains(xRef,'.nii'); xRef=xRef(1:end-4);end
        [x,~,MT]=readNIIExt(xRef, {''});%Read NII
        x=x{1};MT=MT{1}; xRef=[];%Assign array and s_form rotation matrix
    else; x = xRef;xRef=[];
    end
end

%% Make parameters exhaustive
N = size(x);N(end+1:12)=1;
if ~isreal(x);x=abs(x); warning('plotND:: Absolute taken from complex data (x).'); end

%%% Slice information
if ~iscell(sliceInfo); sliceInfo = {sliceInfo, 2};end%make slices stacked horizontally
if isempty(sliceInfo{1}); sliceInfo{1}=1:3;end%Plot all 3 slices
assert( ~any(~ismember( sliceInfo{1} ,[1,2,3])),'plotND:: Invalid slice display (must be sagital (sliceInfo{1}=1), coronal (sliceInfo{1}=2) and/or transversal (sliceInfo{2}=3).');

if length(sliceInfo)<2 || isempty(sliceInfo{2}); sliceInfo{2}=2;end%Plot slices horizontally
assert( ismember( sliceInfo{2} ,[1,2]),'plotND:: Invalid stacking of slices (must be vertically (sliceInfo{2}=1) or horizontally (sliceInfo{2}=2).');

%%% limPlot
if isempty(limPlot);limPlot=[min(x(:))  max(x(:))];end
if length(limPlot)>2; limPLotRef = limPlot(3:4);limPlot = limPlot(1:2);else limPLotRef=[];end%Limits for reference contained in elemts 3:4

%%% Colormap
if isscalar(cmapInfo)
    if cmapInfo==0;cMap = gray(200);
    elseif cmapInfo== 1; cMap = jet(200);
    elseif cmapInfo== 2; cMap = turbo(200);%Introduced in 2020
    elseif cmapInfo== 3; cMap = hot(200);
    elseif cmapInfo== 4; cMap = cmrmap(200);
    else; error('plotND:: cmapInfor scalar not implemented.'); 
    end
end
if ismatrix(cmapInfo) && size(cmapInfo,1)>1; cMap = cmapInfo; assert(size(cMap,2)==3,'plotND:: cmapInfo matrix containing colormap not valid.');end
if ischar(cmapInfo); cMap = colormap(cmapInfo);end
if iscell(cmapInfo) ; cMap = createColormap(cmapInfo);end
useColor =  ~ ( isscalar(cmapInfo) && cmapInfo==0);

%%% Text
if isempty(labelText);labelText=cell(N(4),N(5));end

%%% Set slice indices
if ~iscell(idPlot)
    isRASId=0;
else
    if length(idPlot)<2; isRASId=0;else isRASId=idPlot{2};end
    idPlot = idPlot{1};
end

if isempty(idPlot);idPlot=ceil((N+1)/2);idPlot(4:end)=[];end
idPlot = round(idPlot);%Round for safety
if length(idPlot)==1
    idTemp = ceil((N+1)/2); idTemp(4:end)=[];
    idTemp(1)=idPlot; idPlot=idTemp;idTemp=[];
end

%% Add a referecen image (xRef) - and stack it in the right dimension
colorRef = 1;
if ~isempty(xRef)
    if iscell(xRef)
        if length(xRef)>1 && ~isempty(xRef{2});dimPlotRef = xRef{2};else; dimPlotRef=1; end 
        if length(xRef)>2 && ~isempty(xRef{3});colorRef = xRef{3};else; colorRef=0; end
        xRef=xRef{1};
    else
        dimPlotRef=1;colorRef = 0;%By default, store in the horizontal way and keep the grey colorbar
    end
    if~isreal(xRef);xRef=abs(xRef);warning('plotND:: Absolute taken from complex reference data (xRef).');end%%Take magnitude
    NxRef = size(xRef); NxRef(end+1:12)=1;
    if dimPlotRef==1; dimStoreRef = 4; dimRepRef=5; elseif  dimPlotRef==2; dimStoreRef = 5;dimRepRef=4; end 
    if NxRef(dimRepRef)~=N(dimRepRef); Nrepmat = ones(1,12); Nrepmat(dimRepRef) = N(dimRepRef);xRef=repmat(xRef,Nrepmat);end
    xRef = rescaleND(xRef,limPlot,limPLotRef);
    x = cat(dimStoreRef,xRef,x);
    if size(labelText,mod(dimStoreRef,2)+1) < size(x,dimStoreRef)
        if dimStoreRef==4
            labelText = cat(mod(dimStoreRef,2)+1, cell(1,N(5)), labelText);
        else
            labelText = cat(mod(dimStoreRef,2)+1, cell(N(4),1), labelText);
        end
    end
end

%% Permute to have in RAS (Right-Anterior-Superior) to have consistent imaging planes
if ~isempty(MT)
    [perm, fl] = T2perm(MT(1:3,1:3));%Flips and permutes applied to the array to make logical axes approximate the RAS axes.
    x = flipPermute(x,fl,perm,1); %forward to go to RAS
    if usageM{1}>0 && ~isempty(M); M = flipPermute(M,fl,perm,1);end
    if isRASId==0
        idPlot(fl==1) = (N(fl==1)+1) - idPlot(fl==1);    
        idPlot = idPlot(perm(1:3));
    end
end

%%% Update sizes of the appropriate arrays/structures
N = size(x);
N(end+1:12) = 1;
Nt = size(labelText); 
labelText = dynInd(labelText, {(Nt(1)+1):N(4), (Nt(2)+1):N(5)},1:2, {[]});

%% Stack image planes together
Ni = max(N(1:3)) * ones(1,2);%Don't resample but use margins (avoids Gibsringing artefacts)
X = []; %Array with the stacked images
Xm=[]; %Array with the stacked masks
XmGrey = []; %Array with the stacked masks (the ones that should be gryy for the reference )

for vert = 1:N(4)
    r = [];rm = [];rmGrey = [];
    for hor=1:N(5)
        w=[];wm = [];
        y = dynInd(x,[vert hor],4:5);
        if usageM{1}==1 && ~isempty(M); y = M.*y;end%Multiplicative mask
        
        for d=sliceInfo{1} %which slices to plot
            z=shiftdim(y,d);
            if usageM{1}>=2; m=shiftdim(M,d);end%For contour plotting or color masking
            
            idShift = circshift(idPlot,-d);
            z=dynInd(z, idShift(3)  ,3);
            if usageM{1}>=2;m=dynInd(m, idShift(3)  ,3);end
            
            %%% Apply shifting to get right view - assuming x has approximately LR-PA-IS as dimensions
            if d==1; z = z.';z = flip(z,1); if usageM{1}>=2; m = m.';m = flip(m,1); end ; end
            if d==2; z = flip(z,1); if usageM{1}>=2;m = flip(m,1); end ; end
            if d==3; z = z.'; z = flip(z,1); if usageM{1}>=2; m = m.'; m = flip(m,1); end ; end
            
            z = resampling(z,Ni,2);    
            if usageM{1}>=2; m = single ( resampling(m,Ni,2) > 0.5);end
            
            w=cat(sliceInfo{2},w,(z)); %in which direction to stack different slices
            if usageM{1}>=2; wm = cat(sliceInfo{2}, wm, m);end
        end
        r=cat(2,r,w);
        if usageM{1}>=2 ; rm = cat(2, rm, wm);end
        if ~colorRef && ((dimStoreRef==4 && vert==1)||(dimStoreRef==5 && hor==1)) ; rmGrey = cat(2, rmGrey, ones(size(w),'like',real(w))); else; rmGrey = cat(2, rmGrey, zeros(size(w),'like',real(w)));end %w here for case useM=0 and colorReg=0
    end
    X = cat(1,X,r);
    if usageM{1}>=2 ; Xm = cat(1, Xm, rm);end
    if ~colorRef ; XmGrey = cat(1, XmGrey, rmGrey);end
end

%% Plot
if ~(nargout > 2) %If stacked array with images is requested, don't plot
    
    %%% FIGURE
    if ~isempty(numFig);h=figure(numFig);clf;else; h = figure();end   %set(axes,'Position',[0.005,0.05,0.95,.95]);
    
    %%% IMAGE
    if useColor && usageM{1}==2
        plotMasked(X, limPlot, cMap, [], XmGrey,[], usageM{1}==4);       
    elseif useColor
        plotMasked(X, limPlot, cMap, Xm, XmGrey,[], usageM{1}==4);
    elseif usageM{1}==5 %Transparant mask
        if length(usageM)<2 || isempty(usageM{2}); usageM{2} = .5;end%Opacity
        rgbMask = cat(3,onesL(X),(1-usageM{2}.*Xm),(1-usageM{2}.*Xm));
        X = bsxfun(@times, rescaleND(X,[0 1]), rgbMask);
        imshow(X,[]); if~isequal( limPlot,[0 0]);caxis(gather(limPlot));end
    else
        imshow(X,limPlot);
    end
    colorbar; colormap(cMap);%If you do this after the text, it won't be applied
    set(h,'color','w','Position',get(0,'ScreenSize'));
    
    %%% LABEL TEXT
    textOffset = [10,10];%In x and y axes (horizontal and vertical) - MATLAB convention - this might not generalise well if array sizes are very different
    labelFontSize = FontSizeScaling(2)*25; Linewidth = 15;
    titleFontSize = FontSizeScaling(1)*20;
    colorbarFontSize = FontSizeScaling(3)*20;
    if ~useColor || usageM{1}==3; color='w'; else; color = [0 0 0]; end %black or white
    for i=1:N(4)
        for j = 1:N(5)
            if isempty(labelText{i,j});labelText{i,j}=''; end 
            text(textOffset(1) + (j-1) * Ni(2) * (1+ (sliceInfo{2}==2)*(length(sliceInfo{1})-1)) ,...
                 textOffset(2) + (i-1) * Ni(1) * (1+ (sliceInfo{2}==1)*(length(sliceInfo{1})-1)) ,...
                 labelText{i,j},...
                 'FontSize',labelFontSize,'Linewidth',Linewidth,'Color',color,'Interpreter','latex')
        end
    end
    %hText = findobj(gcf,'Type','Text');
    
    %%% COLORBAR
    hColorbar = findobj(gcf,'Type','Colorbar'); 
    if ~isempty(colorbarText);set(get(hColorbar,'label'),'string',colorbarText,'FontSize',colorbarFontSize,'Interpreter','latex','Linewidth',Linewidth);end
    set(hColorbar,'TickLabelInterpreter','latex','FontSize',colorbarFontSize);
    
    %%% CONTOUR
    if usageM{1}==2
        if length(usageM)<2 || isempty(usageM{2}); c = 'r';else c = usageM{2};end%Colour of contour
        if length(usageM)<3 || isempty(usageM{3}); Linewidth = 2;else Linewidth = usageM{3};end%Linewidth of contour
        if length(usageM)<4 || isempty(usageM{4}); sp = 0.2;else sp = usageM{4};end %For smoothing contour (from 0 to 1 with 1 almost no effect and 0 fully smoothed)
        
        HY = buildFilter(2*size(Xm),'tukeyIso',sp,0,0.3,1);%In DCT domain
        Xm = single( filtering(Xm, HY, 1) > 0.5);%Make contour smoother
        B = bwboundaries(gather(Xm==1));
        
        figure(h); hold on
        for k = 1 : length(B)
             thisBoundary = B{k};  % Get k'th boundary
             x = thisBoundary(:, 2); y = thisBoundary(:, 1);
             plot(x, y, 'color',c, 'LineWidth', Linewidth);
        end
        hold off;
    else
        B=[];
    end
    X=[]; Xm=[];
    
    %%% TITLE
    if ~isempty(titleText);title(titleText,'FontSize',titleFontSize,'Interpreter','latex');end
else
    h=[];hColorbar=[];
end

%%% SAVE
if ~isempty(folderName) && ~isempty(fileName)
    if ~exist(folderName,'dir');mkdir(folderName);end
    saveFig(strcat(folderName,filesep,fileName), [], [], [], 0, [], [], '.png');
end

