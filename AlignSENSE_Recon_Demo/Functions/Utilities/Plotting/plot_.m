function [] = plot_(x, rrx, id, f, fInt,rrf, rrError, usejet, M, useM, numFig)

N = size(x);

if nargin<2 || isempty(rrx); rrx = [];end
if nargin<3 || isempty(id); id = ceil((N+1)/2);end
if nargin<4 || isempty(f); f = [];end
if nargin<5 || isempty(fInt); fInt = [];end
if nargin<6 || isempty(rrf); rrf = [];end
if nargin<7 || isempty(rrError); rrError = [];end
if nargin<8 || isempty(usejet); usejet = 0;end
if nargin<9 || isempty(M); M = [];end
if nargin<10 || isempty(useM); useM=0;end
if nargin<11 || isempty(numFig); numFig=222;end
%%% set slice index
if length(id)==1
    idtemp = ceil((N+1)/2); 
    idtemp(1)=id;id=idtemp;
end
%%% apply masking
if ~isempty(M) && useM
    fInt = fInt.*M;
    f= M.*f;
end

if ~isempty(f) && ~isempty(fInt); error = abs(f-fInt);else error =[];end
nRows = ~isempty(x) + ~isempty(f) + ~isempty(fInt)+ ~isempty(error);

%%% Plot
h=figure(numFig);
set(h,'color','w');

subplot(nRows,3,1); slice = squeeze(x(id(1),:,:));  imshow((slice),rrx);
    colorbar ; cmap = gray(1000); colormap(gca, cmap)
subplot(nRows,3,2); slice = squeeze(x(:, id(2),:));  imshow((slice),rrx);
    colorbar ;cmap = gray(1000); colormap(gca, cmap)
subplot(nRows,3,3); slice = squeeze(x(:,:,id(3)));   imshow((slice),rrx);
    colorbar ;cmap = gray(1000); colormap(gca, cmap)
    
if ~isempty(f)
    subplot(nRows,3,4); slice = squeeze(f(id(1),:,:));  imshow(slice,rrf);
        colorbar; title(''); if usejet; colormap(jet(1000));end
    subplot(nRows,3,5); slice = squeeze(f(:, id(2),:));  imshow(slice,rrf);
        colorbar ; if usejet;colormap(jet(1000));end
    subplot(nRows,3,6); slice = squeeze(f(:,:,id(3)));   imshow(slice,rrf);
        colorbar ; if usejet; colormap(jet(1000));end
    if ~isempty(fInt)
        subplot(nRows,3,7); slice = squeeze(fInt(id(1),:,:));  imshow(slice,rrf);
            colorbar; title('') ;if usejet; colormap(jet(1000));end
        subplot(nRows,3,8); slice = squeeze(fInt(:, id(2),:));  imshow(slice,rrf);
            colorbar ; if usejet; colormap(jet(1000));end
        subplot(nRows,3,9); slice = squeeze(fInt(:,:,id(3)));   imshow(slice,rrf);
            colorbar ; if usejet; colormap(jet(1000));end

        if ~isempty(error)
            subplot(nRows,3,10); slice = squeeze(dynInd( abs(f - fInt),id(1),1));   imshow(slice,rrError);
                colorbar; title('Error ') ; if usejet; colormap(jet(1000));end
            subplot(nRows,3,11); slice = squeeze(dynInd( abs(f - fInt),id(2),2));   imshow(slice,rrError);
                colorbar ; if usejet; colormap(jet(1000));end
            subplot(nRows,3,12); slice = squeeze(dynInd( abs(f - fInt),id(3),3));   imshow(slice,rrError);
                colorbar ; if usejet; colormap(jet(1000));end
        end
    end
end

