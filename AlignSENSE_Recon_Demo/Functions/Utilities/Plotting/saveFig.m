
function [] = saveFig(nameFig, pauseFlag, exportFlag, pauTime, makeTransparant, figResol, suppressWarning, figType) 

%SAVEFIG is a wrapper function to save a generated figure with additional
%options. Code relies on https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig.
%
%   [] = SAVEFIG({NAMEFIG},{PAUSEFLAG},{EXPORTFLAG},{PAUTIME},{MAKETRANSPARANT},{FIGRESOL},{SUPPRESSWARNING},{FIGTYPE})
%   * {NAMEFIG} is the name of the filename to save.
%   * {PAUSEFLAG} is a flag to pause so that user can make changes to the figure.
%   * {EXPORTFLAG} is a flag to block figure exporting.
%   * {PAUTIME} is the time to pause when PAUSEFLAG==1.
%   * {MAKETRANSPARANT} is a flag to make the figure background transparant.
%   * {FIGRESOL} is the resolution of the figure to save.
%   * {SUPPRESSWARNING} is a flag to suppress warnings when calling the export_figure.m function.
%   * {FIGTYPE} is the extension of the figure to save (e.g. png or tiff).
%  
%   Yannick Brackenier 2022-01-30

if nargin < 1 || isempty(nameFig); exportFlag = 0;  end %Disable exporting
if nargin < 2 || isempty(pauseFlag); pauseFlag=0;end
if nargin < 3 || isempty(exportFlag); exportFlag=1;end
if nargin < 4 || isempty(pauTime); pauTime=0;end
if nargin < 5 || isempty(makeTransparant); makeTransparant=0;end
if nargin < 6 || isempty(figResol); figResol=100;end
if nargin < 7 || isempty(suppressWarning); suppressWarning=1;end
if nargin < 8 || isempty(figType); figType= '.png';end %Use .tiff for publications
makeNative = 0;

if suppressWarning; warning('off');end
if exportFlag ==1  
    % Pause - let user change figure size and make changes
    if pauseFlag && pauTime>0; pause(pauTime);elseif pauseFlag; pause();end 
    
    % File names and directory check
    if ~isequal(nameFig(end-4:end), figType); nameFig = strcat(nameFig, figType);end %Check for suffix
    if ~isfolder( fileparts(nameFig)); mkdir(fileparts(nameFig)) ;end%Check if folder exists
    
    % Create input cell array
    varargin={nameFig};%Set filename
    varargin = cat(2,varargin,{sprintf('-r %d',figResol)});
    if makeTransparant; varargin = cat(2,varargin,{'-transparent'}); end
    if makeNative; varargin = cat(2,varargin,{'-native'}); end
    %-------Add extra options here----------
    
    % Export
    export_fig(varargin{:});
    
end

if suppressWarning; warning('on');end

