function files = saveFigs(varargin)
%SAVEFIGS Save open figures to files.
%   SAVEFIGS() saves all the opened figures into the folder ".\fig" in PNG
%   format.
%
%   F = SAVEFIGS() returns a cell array containing the list of the created
%   files (in relative path).
%   
%   SAVEFIGS(DIR) saves  all the opened figures into the folder specified 
%   in DIR as relative path in PNG format.
%
%   SAVEFIGS(DIR,options) saves  all the opened figures into the folder 
%   specified in DIR as relative path with options specified as Param/Value
%   pairs.
%
%   NOTE: I've changed DIR to operate as absolute path
%
%   Options:
%   - 'format'  export format, between PNG (default), JPG, FIG, PDF.
%   - 'handle'  array or single handle to figures to save. Without this
%   argument, all the open figures are saved.
%   - 'style'   name of the custom export style, created in
%       Figure:File->Export Setup.
%   - 'name'    name given to the file if figure.Name property is not present or
%   invalid. Default is 'fig'.
%   - 'dpi'     figure resolution in dots per inch (DPI). Default is 150
%   
%   Notes: If the default or specified folder is not present, it will be
%   created. In case a file with the same name is present in the folder,
%   the new file will be created with a trailing progressive number, in
%   order not to overwrite it.

%   Copyright (c) 2015, Edoardo Bezzeccheri
%   All rights reserved.

% Initialize the input argument parser
p = inputParser;
defaultDir = 'img';
defaultFormat = 'png';
expectedFormats = {'png','fig','jpg','pdf'};%TODO: add more formats
defaultStyle = '';
defaultFigureName = 'fig';
defaultPdfDpi = 150;
defaultHandle = '';

addOptional(p,'dir',defaultDir,@isstr);
addParameter(p,'format',defaultFormat,...
    @(x) any(validatestring(x,expectedFormats)));
addParameter(p,'handle',defaultHandle,...
    @(x) any(ishandle(x)));
addParameter(p,'style',defaultStyle,@isstr);
addParameter(p,'name',defaultFigureName,@isstr);
addParameter(p,'dpi',defaultPdfDpi,@isnumeric);

parse(p,varargin{:});

% Check for open figures
if(~isempty(p.Results.handle))
    figures = p.Results.handle;
else
    root = groot;
    figures = root.Children;
end
assert(~isempty(figures),'No open figures.');

figdir = p.Results.dir;
if(isdir(figdir))
    disp('Folder already present, saving there');
else
    mkdir(figdir);
    disp(['Created folder ' figdir]);
end

isnamevalid = @(x) isempty(regexp(x, '[/\*:?"<>|]', 'once'));
for i=1:length(figures);
    fig = figures(i);
    
    % File name creation
    figname = fig.Name;
    if(~isnamevalid(figname) || isempty(figname) || length(figname) > 60)
        figname = [p.Results.name num2str(fig.Number)];
        warning('Figure (%d) has empty or invalid name for file creation. Using its number instead.',...
            fig.Number);
    end
    %filepath = fullfile('.',figdir,figname);  
    filepath = fullfile(figdir,figname);  
    
    % Changing style (if needed)
    if (~isempty(p.Results.style))
        style = hgexport('readstyle',p.Results.style);
        hgexport(figures(i),'temp_dummy',style,'applystyle', true);
    end
    
    % Saving file
    switch p.Results.format
    case 'png'
        filepath = checkFilename([filepath '.png']);
        print(fig,'-dpng',filepath);
    case 'fig'
        filepath = checkFilename([filepath '.fig']);
        savefig(fig,filepath);
    case 'jpg'
        filepath = checkFilename([filepath '.jpg']);
        print(fig,'-djpeg',filepath);
    case 'pdf'
        filepath = checkFilename([filepath '.pdf']);
        printPdf(fig,filepath,p.Results.dpi);
    end
    
    % Reverting to factory style (only if a style was indicated)
    if (~isempty(p.Results.style))
        style = hgexport('factorystyle');
        hgexport(figures(i),'temp_dummy',style,'applystyle', true);
    end
    
    % Returning paths
    files{i} = filepath;
end

end

% Recursive function needed for checking if there is another file with the
% same name. In that case it will add a number at the end.
function [filepath] = checkFilename(filepath,varargin)

    [dir,filename,ext] = fileparts(filepath);
    
    if (nargin == 2)
        rec = varargin{1};
        
        k = strfind(filename,'-');
        figname = filename(1:k(end)-1);
    else 
        rec = 1;
        figname = filename;
    end
    
    filepath = [fullfile(dir,filename) ext];
    if(exist(filepath,'file') == 2)
        filename = [figname '-' num2str(rec)];
        filepath = checkFilename([fullfile(dir,filename) ext],rec+1);
        
        if(rec == 1) 
            warning('File with same (%s) name present in the chosen directory. Using unique name.',...
        [figname ext]);
        end
    end
    
end

% The code below mostly belongs to Gabe Hoffmann [gabe.hoffmann@gmail.com],
% which has written the function "save2pdf", originally available at 
% http://www.mathworks.com/matlabcentral/fileexchange/16179-save2pdf
function printPdf(fig,filepath,dpi)
    % Backup previous settings
    settings = {'PaperType','PaperUnits','Units',...
        'PaperPosition','PaperSize','PaperOrientation'};
    preSettings = get(fig,settings);
    
    % Make changing paper type possible
    set(fig,'PaperType','<custom>');
    set(fig,'PaperOrientation','landscape');
    
    % Set units to all be the same
    set(fig,'PaperUnits','inches');
    set(fig,'Units','inches');

    % Set the page size and position to match the figure's dimensions
    position = get(fig,'Position');
    set(fig,'PaperPosition',[0,0,position(3:4)]);
    set(fig,'PaperSize',position(3:4));
    
    % Save the pdf (this is the same method used by "saveas")
    print(fig,'-dpdf',filepath,sprintf('-r%d',dpi))
    
    % Restore the previous settings
    set(fig,settings,preSettings);
end
