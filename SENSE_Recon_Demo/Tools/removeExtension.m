function fil=removeExtension(fil,ext)

%REMOVEEXTENSION   Removes a series of extensions from a file name (or a
%cell array of filenames)
%   FIL=REMOVEEXTENSION(FIL,{EXT})
%   * FIL is the file name whose extensions we want to remove
%   * {EXT} is a cell (or plain char vector) of extensions that this file 
%   name should not have
%   ** FIL is the resulting file name
%

if ~iscell(fil)  
    while 1
        [path,file,exte]=fileparts(fil);      
        if any(strcmp(exte,ext));fil=fullfile(path,file);else break;end
    end
else
    for n=1:length(fil);fil{n}=removeExtension(fil{n},ext);end
end

