

function [xRGB] = plotMasked(x, rr, cmap, M, Mgrey, NColors, makeTransparant)

%PLOTMASKED converts a greyscale array into a RGB array corresponding to a colormap but with a black background specified by a mask.
%
%   [XRGB] = PLOTMASKED(X,{RR},{CMAP},{M},{MGREY},{NCOLORS},{MAKETRANSPARANT})
%   * X is the greyscale array.
%   * {RR} is the dynamic range of the greyscale array.
%   * {CMAP} is the colormap to use.
%   * {M} is the mask to use to set a black background.
%   * {MGREY} is the mask to use to keep greyscale values.
%   * {NCOLORS} is the number of colors to use in the colormap.
%   * {MAKETRANSPARANT} whether to make the background transparant instead of black.
%   ** XRGB is the output array in RGB values for the given result that can be called with imshow.m.
%
%   Yannick Brackenier 2022-01-30


if ~isreal(x); warning('plotMasked:: Image should not be complex.\n Color assignments will fail so real part taken!'); x = real(x);end

if nargin < 2 || isempty(rr); rr = gather([min(x(:)), max(x(:))]);end
if nargin < 3 || isempty(cmap); cmap = [];end
if nargin < 4 || isempty(M); M = ones(size(x));end%Mask to set to black ([0 0 0] RGB color)
if nargin < 5 || isempty(Mgrey); Mgrey = [];end%Mask to set which voxels should remain grey in image
if nargin < 6 || isempty(NColors); NColors = 1000;end%number of contrasts (levels) in the colormap
if nargin < 7 || isempty(makeTransparant); makeTransparant = 0;end

%cmap = jet(Ncolors+1);
if isempty(cmap);cmap = jet(NColors+1);else; NColors=size(cmap,1)-1;end 
cmap = permute(cmap,[1 3 2]);

%%% Constrain range of image
x (x<rr(1)) = rr(1);
x (x>rr(2)) = rr(2);

%%% Rescale to number of levels
if ~isempty(Mgrey)
    x(Mgrey==0) = round(rescaleND(x(Mgrey==0),[0,NColors], rr));
    x(Mgrey==1) = rescaleND(x(Mgrey==1),[0 1], rr); %Voxels that should remain gray must be rescaled to [0 1] without rounding
else%all voxels in color
    x = round(rescaleND(x,[0,NColors], rr));
end
Nx = size(x);Nx(end+1:12) = 1;

%%% Make RGB version and reshape
xRGB = repmat(x, [1 1 3]);
NxRGB = size(xRGB);NxRGB(end+1:12) = 1;

x = reshape(x, [prod(Nx(1:2)) 1 Nx(3)]);
xRGB = reshape(xRGB, [prod(NxRGB(1:2)) 1 NxRGB(3)]);
if ~isempty(Mgrey); Mgrey = reshape(Mgrey, [prod(NxRGB(1:2)) 1 ]);end

%%% Assign colors
for i=(0:NColors)
    idx = x==i;
    if ~isempty(Mgrey); idx = (idx==1 & Mgrey==0);end%remove indices of voxels that should remain grey
    xRGB = dynInd( xRGB, idx,1, repmat(dynInd(cmap, i+1, 1),[numel( find(idx)) 1])) ; 
end

%%%Reshape back and set M black ([0 0 0])
xRGB = reshape(xRGB, NxRGB(1:3));
xRGB = bsxfun(@times, xRGB, M);%Applies for all pixels (colours and grey)

%%% Plotting
if nargout < 1
    if makeTransparant
        hIm = imshow(xRGB,[]); 
        set(hIm,'AlphaData',gather(M));
    else
        h=gcf; figure(h);
        imshow(xRGB,[]);
    end
    colorbar ; if~isequal( rr,[0 0]);caxis(gather(rr));end
    xRGB = [];
end

