

function [xRGB] = plotMaskedOld(x, rr, M, Mgrey, Ncolors, makeTransparant)
%%% plot x with range rr in jet colormap and set pixels in M to black

if ~isreal(x); warning('plotMasked:: Image should not be complex.\n Color assignments will fail so real part taken!'); x = real(x);end

if nargin < 2 || isempty(rr); rr = gather([min(x(:)), max(x(:))]);end
if nargin < 3 || isempty(M); M = ones(size(x));end%Mask to set to black ([0 0 0] RGB color)
if nargin < 4 || isempty(Mgrey); Mgrey = [];end%Mask to set which voxels should remain grey in image
if nargin < 5 || isempty(Ncolors); Ncolors = 1000;end%number of contrasts (levels) in the jet colormap
if nargin < 6 || isempty(makeTransparant); makeTransparant = 0;end%number of contrasts (levels) in the jet colormap

cmap = jet(Ncolors+1);
cmap = permute(cmap,[1 3 2]);

%%% Constrain range of image
x (x<rr(1)) = rr(1);
x (x>rr(2)) = rr(2);

%%% Rescale to number of levels
if ~isempty(Mgrey)
    x(Mgrey==0) = round(rescaleND(x(Mgrey==0),[0,Ncolors], rr));
    x(Mgrey==1) = rescaleND(x(Mgrey==1),[0 1], rr); %Voxels that should remain gray must be rescaled to [0 1] without rounding
else%all voxels in color
    x = round(rescaleND(x,[0,Ncolors], rr));
end
Nx = size(x);Nx(end+1:12) = 1;

%%% Make RGB version and reshape
xRGB = repmat(x, [1 1 3]);
NxRGB = size(xRGB);NxRGB(end+1:12) = 1;

x = reshape(x, [prod(Nx(1:2)) 1 Nx(3)]);
xRGB = reshape(xRGB, [prod(NxRGB(1:2)) 1 NxRGB(3)]);
if ~isempty(Mgrey); Mgrey = reshape(Mgrey, [prod(NxRGB(1:2)) 1 ]);end

%%% Assign colors
for i=(0:Ncolors)
    idx = x==i;
    if ~isempty(Mgrey); idx = (idx==1 & Mgrey==0);end%remove indexes of voxels that should remain grey
    xRGB = dynInd( xRGB, idx,1, repmat(dynInd(cmap, i+1, 1),[numel( find(idx)) 1])) ; 
end

%%%Reshape back and set M black ([0 0 0])
xRGB = reshape(xRGB, NxRGB(1:3));
xRGB = bsxfun(@times, xRGB, M);%Applies for all pixels (colors and grey)

%%% Plotting
if nargout < 1
    if makeTransparant
        hIm = imshow(xRGB,[]); 
        set(hIm,'AlphaData',M);
    else
        imshow(xRGB,[]);
    end
    colorbar ; if~isequal( rr,[0 0]);caxis(gather(rr));end
    xRGB = [];
end