%% Shortcut examples
% I find very useful to use this function as a toolbar shortcut. By using
% the provided code for the shortcut, it saves all the opened figures in
% the folder specified in the variable "figdir". 
% If the variable is not present, it saves them into the default folder .\fig.

% The example below can be used directly for shortcuts.

%% Basic Shortcut
disp('Saving open figures...');
if(exist('figdir','var') == 1)
    saveFigs(figdir);
else
    saveFigs();
end

%% Shortcut for PDF prints (300 DPI)
disp('Saving open figures into PDFs @ 300 DPI...');
if(exist('figdir','var') == 1)
    saveFigs(figdir,'format','pdf','dpi',300);
else
    saveFigs('img','format','pdf','dpi',300);
end

%% Shortcut for MATLAB figures named 'measurement-#'
disp('Saving open figures into MATLAB .fig...');
if(exist('figdir','var') == 1)
    saveFigs(figdir,'format','fig','name','measurement');
else
    saveFigs('img','format','fig','name','measurement');
end
