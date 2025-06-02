function findString(targString,openFil,path,subdirs)

%FINDSTRING   Looks for a given string in a set of files
%   FINDSTRING(TARGSTRING,{OPENFIL},{PATH},{SUBDIRS})
%   * TARGSTRING is the string to be looked for
%   * {OPENFIL} is a flag to open the files where occurrences have been 
%   found in a matlab editor. Usefool when looking for code. Defaults to 0
%   * {PATH} is the path (directory or file) where to look for the string.
%   Defaults to the current matlab path
%   * {SUBDIRS} is a flag to indicate whether to recursively look in
%   subfolders when looking in directories. It defaults to 1
%

if ~exist('openFil','var') || isempty(openFil);openFil=0;end
if ~exist('path','var') || isempty(path);path=pwd;end
if ~exist('subdirs','var') || isempty(subdirs);subdirs=1;end

structPath = dir(path);%Returns all the files and folders in the directory
structPath(ismember( {structPath.name}, {'.', '..'}))=[];

if subdirs
    structSubDir=structPath([structPath.isdir]);
    for n=1:length(structSubDir);findString(targString,openFil,fullfile(path,structSubDir(n).name),subdirs);end%Recursive call to subdirs
end

structFil=structPath(~[structPath.isdir]);
for n=1:length(structFil)
    if ~strcmp(structFil(n).name(end),'~')
        fileName=fullfile(path,structFil(n).name);
        foundInFile=0;
        fid=fopen(fileName);
        m=1;
        while ~feof(fid)
           stringLine=fgetl(fid);
           found=strfind(stringLine,targString);  
           if ~isempty(found);fprintf('File %s. Line %d. Text: "%s"\n',fileName,m,stringLine);foundInFile=1;end
           m=m+1;
        end
        fclose(fid);       
        if openFil && foundInFile;open(fileName);end%To inspect the results in a matlab editor when looking for code
    end
end