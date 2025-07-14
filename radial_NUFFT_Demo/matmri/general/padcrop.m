function y = padcrop(x,szIn,dim)
%
%  res = padcrop(x,szIn,dim)
%  Changes size of x along dimension "dim" by either zero padding or cropping around its center
%
% (c) Corey Baron 2015
% 

if numel(szIn)>1 && nargin>2
    error('Incompatible inputs: with dim supplied, szIn should be an integer')
end
if nargin<3
    dim = [];
end

if ~isempty(dim)
    y = singDim(x,szIn,dim);
else
    for dim=1:length(szIn)
        x = singDim(x,szIn(dim),dim);
    end
    y = x;
end

end
    

function y = singDim(x,szIn,dim) 
% Determine number to pad/crop
sz = [size(x), ones(1,max(ndims(x),dim))];
npc = szIn-sz(dim);

if npc > 0
    % Set up inds for defining size of zeros array
    inds = true(size(sz));
    inds(dim:end) = false;
    inds2 = ~inds;
    inds2(dim) = false;
    % Finish padding
    if (dim <= 3) && (ndims(x) < 11)
        % Here we do a less computationally expensive approach for the most
        % common cases
        szx = [size(x),ones(1,5)];
        if dim==1
            y = zeros([szIn szx(2:end)], 'like', x);
            y(floor(npc/2)+1:floor(npc/2)+szx(dim),:,:,:,:,:,:,:,:,:) = x;
        elseif dim==2
            y = zeros([szx(1) szIn szx(3:end)], 'like', x);
            y(:,floor(npc/2)+1:floor(npc/2)+szx(dim),:,:,:,:,:,:,:,:) = x;
        elseif dim==3
            y = zeros([szx(1:2) szIn szx(4:end)], 'like', x);
            y(:,:,floor(npc/2)+1:floor(npc/2)+szx(dim),:,:,:,:,:,:,:) = x;
        end
    else
        zpadding1 = zeros([sz(inds),floor(npc/2),sz(inds2)], 'like', x);
        zpadding2 = zeros([sz(inds),npc-floor(npc/2),sz(inds2)], 'like', x);
        y = cat(dim,zpadding1,x,zpadding2);
    end
elseif npc<0
    npc = -npc;
    % Finish cropping
    if (dim <= 3) && (ndims(x) < 11)
        % Here we do a less computationally expensive approach for the most
        % common cases
        if dim==1
            y = x(floor(npc/2)+1:floor(npc/2)+szIn,:,:,:,:,:,:,:,:,:);
        elseif dim==2
            y = x(:,floor(npc/2)+1:floor(npc/2)+szIn,:,:,:,:,:,:,:,:);
        elseif dim==3
            y = x(:,:,floor(npc/2)+1:floor(npc/2)+szIn,:,:,:,:,:,:,:);
        end
    else
        % Here we accomodate any other case in a slightly more expensive
        % way
        x = permute(x, [dim 1:dim-1 dim+1:ndims(x)+1]);
        sz_a = size(x);
        y = zeros([szIn, sz_a(2:end)], 'like', x);
        y(:,:) = x(floor(npc/2)+1:floor(npc/2)+szIn,:);
        y = permute(y, [2:dim 1 dim+1:ndims(x)+1]);
    end
else
    y = x;
end
end









