
function [ h, hColorbar,X, B] = plotND(xpre, x, rrx, id, useJet, dList, Text, M, useM, numFig,MT)

%PLOTND is a wrapper function to facilitate plotting multiple volumes at a single comment with additional functionality 
%   (e.g. plot 5D arrays, arrange views flexibly, add text, mask, orient images, etc.).
%   [H, HCOLORBAR,X,B]=PLOTND(XPRE,X,RRX,ID,USEJET,DLIST,TEXT,M,USEM,NUMFIG,MT)
%   * {XPRE} is a cell structure containing information about a referene volume to add. XPRE{1} contains the volume, 
%     XPRE{2} the dimension along which to add the reference and XPRE{3}  whether to use the jet colorbar.
%     If only XPRE is provided, it can consist of a string referrng to a NIFTI file, which will then be read and plotted.
%   * X a 5D array of volumes to plot. The 4th dimension will be plot from left to right and the 5th dimension from top to bottom.
%   * {RRX} the range to plot the arrays. RRX(1:2) the rage for the X array and RRX(3:4) can optionally contain the XPRE range.
%   * {ID} the id of the slices to plot. Defaults to the centre. ID as in original array provided. Can be taken directly from NIFTI viewers (e.g.ITK-SNAP)
%   * {USEJET} flag to use the jet colormap for plotting the array X. Defaults to 0.
%   * {DLIST} a cell array containing information about how to plot the slices. DLIST{1} contains the slices to plot (sag,cor or transverse) and DLIST{2} along which dimension to plot them (defaults to 2=horizontal).
%   * {TEXT} a cell array containing string for the titles used for the subfigures. Size of the array must be consistent with the dimensions of X.
%   * {M} the mask to use for plotting.
%   * {USEM} flag to indicate how to use the mask. 1 if multiplicative, 2 if for contours and 3 if as a black mask for the jet coloured array X.
%   * {NUMFIG} the number of the figure.
%   * {MT} the orientation information to have correct orientation for the plotted slices.
%   ** H the handle of the figure.
%   ** HCOLORBAR the handle of the colorbar.
%   ** X the 2D array containing all the concatenated sub-figures of the slices.
%   ** B the boundary information for when the mask is used for contour plots.
%

%%% Check if one argument provided: can be nifti file or the actial image 
useNIFTI=0;
if nargin==1 || isempty(x)
    if ischar(xpre)
        [x,~,MT]=readNIIExt(xpre, {''});
        x=x{1};MT=MT{1}; xpre=[];
        if~isreal(x);x=abs(x);warning('plotND:: Absolute taken from complex data.');end
        useNIFTI=1;
    else; x = xpre;xpre=[];
    end
end

N = size(x);N(end+1:12)=1;
if ~isreal(x);x=abs(x);warning('plotND:: Absolute taken from complex data.');end
if nargin<3 || isempty(rrx); rrx = [min(x(:))  max(x(:))];end
if nargin<4 || isempty(id); id = ceil((N+1)/2);id(4:end)=[];end
if nargin<5 || isempty(useJet); useJet = 0;end
if nargin<6 || isempty(dList); dList = {1:3,2};end
if nargin<7 || isempty(Text); Text=cell(N(5),N(4));end%5th dimension stored vertically (top to bottom)
if nargin<8 || isempty(M); M = [];end %binary mask to set background black
if nargin<9 || isempty(useM); useM={0};end %0 if not used,1 if applied to image, 2 if used for contours, 3 if used for setting image black (useful if jet color), 4 for transparancy
if nargin<10 || isempty(numFig); numFig=[];end
if (nargin<11 && ~useNIFTI) || isempty(MT); MT=[];end%s_form rotation matrix from NIFTI header (in RAS frame)
colorRef = 1;

if ~iscell(dList); dList = {dList, 1};end %make slices stacked vertically
if isempty(dList{1}); dList{1}=1:3;end ;if isempty(dList{2}); dList{2}=1;end
if length(rrx)>2; rrpre = rrx(3:4);rrx = rrx(1:2);else rrpre=[];end

%%% set slice index
id = round(id);
if length(id)==1
    idtemp = ceil((N+1)/2); idtemp(4:end)=[];
    idtemp(1)=id; id=idtemp;idtemp=[];
end

%%% add a referecen image (xpre) - and stack it in the right dimension
if ~isempty(xpre)
    if iscell(xpre)
        if length(xpre)>1 && ~isempty(xpre{2});dimPlot = xpre{2};else; dimPlot=2; end 
        if length(xpre)>2 && ~isempty(xpre{3});colorRef = xpre{3};else; colorRef=3; end
        xpre=xpre{1};
    else dimPlot=2;colorRef = 0; 
    end%By default, store in the horizontal way and remain grey colorbar
    if~isreal(xpre);xpre=abs(xpre);warning('plotND:: Absolute taken from complex data.');end%%Take magnitude
    Nx1 = size(xpre); Nx1(end+1:12)=1;
    if dimPlot ==1; dimStore = 5;dimRep=4; elseif  dimPlot ==2; dimStore = 4;dimRep=5; end 
    if Nx1(dimRep)~=N(dimRep); Nrepmat = ones(1,12); Nrepmat(dimRep) = N(dimRep);xpre =repmat(xpre,Nrepmat);end
    xpre = rescaleND(xpre,rrx,rrpre);
    x = cat(dimStore,xpre,x);
    if size(Text,mod(dList{2},2)+1) < size(x,dimStore) ; Text = cat(mod(dList{2},2)+1, cell(1+(N(5)-1)*(dList{2}==1),1+(N(4)-1)*(dList{2}==2)), Text);end
end

%%% Permute to have in RAS
if ~isempty(MT)
    [perm, fl] = T2perm(MT(1:3,1:3));%Flips and permutes applies to array compares to RAS reference
    x = flipPermute(x,fl,perm,1); %forward to go to RAS
    if useM{1}>0 && ~isempty(M); M = flipPermute(M,fl,perm,1);end
    id(fl==1) = (N(fl==1)+1) - id(fl==1);
    id = id(perm(1:3));%;%id(perm(1:3))= id(1:3);
end
N = size(x);
N (end+1:12) = 1;
Nt = size(Text); 
Text = dynInd(Text, {(Nt(1)+1):N(5), (Nt(2)+1):N(4)},1:2, {[]});

%%% Stack images together
Ni = max(N(1:3)) * ones(1,2);%Don't resample but use margins (avoids Gibsringing artefacts)Ni = [150 150];%size of the images to be stacked
X = [];Xm=[];XmGrey = [];
for vert = 1:N(5)
    r = [];rm = [];rmGrey = [];
    for hor=1:N(4)
        w=[];wm = [];
        y = dynInd(x,[hor vert],4:5);
        if useM{1}==1 && ~isempty(M); y = M.*y;end
        
        for d=dList{1} %which slices to plot
            z=shiftdim(y,d);
            if useM{1}>=2; m=shiftdim(M,d);end%for contour plotting
            
            idShift = circshift(id,-d);
            z=dynInd(z, idShift(3)  ,3);
            if useM{1}>=2;m=dynInd(m, idShift(3)  ,3);end
            
            %%% Apply shifting to get right view - assuming x has LR-PA-IS as dimensions
            if d==1; z = z.';z = flip(z,1); if useM{1}>=2; m = m.';m = flip(m,1); end ; end
            if d==2; z = flip(z,1); if useM{1}>=2;m = flip(m,1); end ; end
            if d==3; z = z.'; z = flip(z,1); if useM{1}>=2; m = m.'; m = flip(m,1); end ; end
            
            z = resampling(z, Ni,2);%z = resampling(z,Ni);            
            if useM{1}>=2; m = single ( resampling(m,Ni,2) > 0.5);end
            
            w=cat(dList{2},w,(z)); %in which direction to stack different slices
            if useM{1}>=2; wm = cat(dList{2}, wm, m);end
        end
        r=cat(2,r,w);
        if useM{1}>=2 ; rm = cat(2, rm, wm);end
        if ~colorRef && ((dimStore==4 && hor==1)||(dimStore==5 && vert==1)) ; rmGrey = cat(2, rmGrey, ones(size(w),'like',real(w))); else; rmGrey = cat(2, rmGrey, zeros(size(w),'like',real(w)));end %w here for case useM=0 and colorReg=0
    end
    X = cat(1,X,r);
    if useM{1}>=2 ; Xm = cat(1, Xm, rm);end
    if ~colorRef ; XmGrey = cat(1, XmGrey, rmGrey);end
end

%%% Plot
if ~ (nargout > 2) %if output is requested, don't plot
    %%%Image
    if ~isempty(numFig);h=figure(numFig);else h = figure();end
    %set(axes,'Position',[0.005,0.05,0.95,.95]);
    
    if ~(useJet==1 && (useM{1}==3||~colorRef))
        if useM{1}==4; hIm = imshow(X,rrx); set(hIm, 'AlphaData',gather(Xm));else; imshow(X,rrx);end
    else
        plotMasked(X, rrx, Xm, XmGrey, 1000, useM{1}==4);
    end
    colorbar; if useJet; cmap = jet(1000); colormap(cmap);end
    set(h,'color','w','Position',get(0,'ScreenSize'));
    
    %%%Text
    textOffset = [10,10];
    FontSize = 20; Linewidth = 15;
    if useJet && useM{1}~=3; color = [0 0 0]; else color='w';end %black or white
    for i=1:N(4)
        for j = 1:N(5)
            if isempty(Text{j,i});Text{j,i}=''; end 
           text(textOffset(1)+ (i-1)*Ni(2),textOffset(2) + (j-1)*Ni(1), Text{j,i} ,...
               'FontSize',FontSize,'Linewidth',Linewidth,'Color',color,'Interpreter','latex')
        end
    end
    %hText = findobj(gcf,'Type','Text');
    
    %%%Colorbar
    hColorbar = findobj(gcf,'Type','Colorbar'); 
    %set(get(hColorbar,'label'),'string','[$Hz/ ^\circ$]','FontSize',FontSize,'Interpreter','latex','Linewidth',Linewidth);
    set(hColorbar,'TickLabelInterpreter','latex','FontSize',FontSize);
    
    %%%Contour
    if useM{1}==2
        if length(useM)<4 || isempty(useM{4}); sp = 0.2;else sp = useM{4};end %For smoothing contour
        HY = buildFilter(2*size(Xm),'tukeyIso',sp,0,0.3,1);
        Xm = single( filtering(Xm, HY, 1) > 0.5);%Make contour smoother
        B = bwboundaries(gather(Xm==1));
        
        figure(h); hold on
        for k = 1 : length(B)
             thisBoundary = B{k};  % Get k'th boundary
             x = thisBoundary(:, 2); y = thisBoundary(:, 1);
             if length(useM)<2 || isempty(useM{2}); c = 'r';else c = useM{2};end
             if length(useM)<3 || isempty(useM{3}); Linewidth = 2;else Linewidth = useM{3};end
             plot(x, y, 'color',c, 'LineWidth', Linewidth);
        end
        hold off;
    else
        B=[];
    end
    X=[]; Xm=[];
else
    h=[];hColorbar=[];
end
