

function cmap = createColormap(colors, m, positions)
  
%CREATECOLORMAP creates a colormap from set of colors that continuously fade
%into each other. Syntax slightly changed from and based on https://uk.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap
%   [CMAP]=CREATECOLORMAP(COLORS, {M},{POSITIONS})
%   * COLORS the colors to be used in HEX, NAME or RGB format. The first
%   color corresponds with the lowest values
%   * {M} the number of samples in the colormap.    
%   * {POSITION} the positions of the colors between 0 and 1.
%   ** CMAP the colormap (dimensions [mx3])
%


if iscell (colors) ; colors = flip( colors );end%YB: function below takes as first argument the highest value
if ~iscell (colors) ; colors = flip( colors,1 );end
% ----------------------------------------------------------------------- %
% FUNCTION "customcolormap" defines a customized colobar given the        %
% positions and the colors that are going to generate the gradients.      %
%                                                                         %
%   Input parameters:                                                     %
%       - positions:    Vector of positions, from 0 to 1. Note that the   %
%                       first position must be 0, and the last one must   %
%                       be 1.                                             %
%       - colors:       Colors to place in each position. This parameter  %
%                       can be specified as a RGB matrix (n_colors x 3), or
%                       as a cell vector, containing HTML values.         %
%                       For instance: {'#ffffff','#ff0000','#000000'} is  %
%                       equivalent to [1 1 1; 1 0 0; 0 0 0].              %
%       - m:            (Optional) Number of points (recommended: m > 64).%
%                                                                         %
%   Output variables:                                                     %
%       - cmap:            Colormap in RGB values (dimensions [mx3]).        %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       cmap = customcolormap([0 0.5 1], {'#ffffff','#ff0000','#000000'});   %
%       colorbar; colormap(cmap);                                            %
% ----------------------------------------------------------------------- %
%   Versions:                                                             %
%       - v1.0.:    (19/11/2018) Original script.                         %
% ----------------------------------------------------------------------- %
%       - Author:   Víctor Martínez-Cagigal                               %
%       - Date:     19/11/2018                                            %
%       - Version:  1.0                                                   %
%       - E-mail:   vicmarcag (at) gmail (dot) com                        %
%                                                                         %
%       Biomedical Engineering Group (University of Valladolid), Spain    %
% ----------------------------------------------------------------------- %
    

    % Error detection and defaults
%     if nargin < 2
%        f = get(groot,'CurrentFigure');
%        if isempty(f)
%           m = size(get(groot,'DefaultFigureColormap'),1);
%        else
%           m = size(f.Colormap,1);
%        end
%     end
    if nargin<2 || isempty(m); m =200;end %YB: take fine enough 
    if ~isnumeric(m), error('Parameter m must be numeric.'); end
    
    if iscell(colors)
        colors = colors(:);
        n_colors = length(colors);
        for i = 1:n_colors
            temp = colors{i};
            if ~ischar(temp)
                error(['Colors must be specified in HEX format (e.g., #FFFFFF).' ...
                ' Type "help colorbar" for further information']);

% YB: disabled since I also allow color names (e.g. 'blue')
%             elseif ~strcmp(temp(1),'#')
%                 error(['Character # is missing if %s.' ...
%                 ' Type "help colorbar" for further information'], temp);
%             elseif length(temp)~=7
%                 error(['Not a valid color format: %s (use this format: #FFFFFF).' ...
%                 ' Type "help colorbar" for further information'], temp);
            end
        end
        hexFormat = strcmp(temp(1),'#');%YB: added
        
    elseif ismatrix(colors)
        n_colors = size(colors);
        if length(n_colors) ~= 2
            error(['RGB colors must be a 2D matrix.' ...
                ' Type "help colorbar" for further information']);
        elseif n_colors(2) ~= 3
            error(['RGB colors matrix must have 3 columns.' ...
                ' Type "help colorbar" for further information']);
        elseif min(colors(:))<0 || max(colors(:))>255
            error(['RGB colors matrix values must range from 0 to 255.' ...
                ' Type "help colorbar" for further information']);
        end
    else
        error(['Colors must be a cell vector or a matrix of RGB values.' ...
            ' Type "help colorbar" for further information']);
    end
    
    if nargin < 3 || isempty(positions); positions = linspace(0,1,n_colors(1)); end%YB

    if ~isvector(positions)
        error(['Positions must be specified as a vector.' ...
            ' Type "help colorbar" for further information']);
    elseif min(positions)<0 || max(positions)>1
        error(['Positions must range from 0 to 1 in an ascending order.' ...
            ' Type "help colorbar" for further information']);
    elseif length(positions) ~= length(unique(positions))
        error(['Check the positions vector, there are some duplicates.' ...
            ' Type "help colorbar" for further information']);
    else
        positions = sort(positions, 'ascend');
        if positions(1)~=0
            error(['The first positions must be 0.' ...
            ' Type "help colorbar" for further information']);
        elseif positions(length(positions))~=1
            error(['The last positions must be 1.' ...
            ' Type "help colorbar" for further information']);
        elseif length(positions) ~= n_colors
            error(['The number of positions does not match the number of colors.' ...
            ' Type "help colorbar" for further information']);
        end
    end
    % Convert HEX OR NAME colors into RGB colors if required
    if iscell(colors)
        if hexFormat
            hex_colors = colors;
            colors = NaN(n_colors,3);
            for i = 1:n_colors
                colors(i,:) = hex2rgb(hex_colors{i});
            end
        else %nameFormat
            name_colors = colors;
            colors = NaN(n_colors,3);
            for i = 1:n_colors
                colors(i,:) = name2rgb(name_colors{i});
            end
        end
    end
    
    % Compute positions along the samples
    color_samples = round((m-1)*positions)+1;
    % Make the gradients among colors
    cmap = zeros(m,3);
    cmap(color_samples,:) = colors;
    diff_samples = diff(color_samples)-1;
    for d = 1:1:length(diff_samples)
        if diff_samples(d)~=0
            color1 = colors(d,:);
            color2 = colors(d+1,:);
            G = zeros(diff_samples(d),3);
            for idx_rgb = 1:3
                g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                g([1, length(g)]) = [];
                G(:,idx_rgb) = g';
            end
            cmap(color_samples(d)+1:color_samples(d+1)-1,:) = G;
        end
    end
    cmap = flipud(cmap);
end
