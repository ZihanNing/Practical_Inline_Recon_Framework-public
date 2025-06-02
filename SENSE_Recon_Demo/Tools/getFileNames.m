function fils=getFileNames(folder,typ)

%GETFILENAMES   Gets a list of filenames on a given folder
%   FILS=GETFILENAMES(FOLDER,{TYP})
%   * FOLDER is the folder name
%   * {TYP} is the type of files. 0 for directories, 1 for files, 2 for 
%   both, it defaults to 1
%   ** FILS is a cell of strings containing the list of files
%

fils=dir(folder);%Returns all the files and folders in the directory
fils(ismember({fils.name},{'.', '..'}))=[];

if typ==1;fils=fils(~[fils.isdir]);elseif typ==0;fils=fils([fils.isdir]);end
fils={fils.name};
