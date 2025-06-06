

function [flag] = existsFileVar(fileName, varName)

%EXISTSFILEVAR checks whether a file and containing variables exist on the file system. 
%   [FLAG] = EXISTSFILEVAR(FILENAME, {VARNAME})
%   * FILENAME is the name of the file of interest.
%   * {VARNAME} is a cell array with the names of the variables of interest.
%   ** FLAG is a list of flags indicating whether the file exists (1) and whether it contains the variables (2).
%
% Yannick Brackenier 2022-07-26

if nargin<1 || isempty(fileName);fileName = '';end
if nargin<2 || isempty(varName);varName = {};end

[~,~,suff] = fileparts(fileName);
if isempty(suff); suff = '.mat';fileName = strcat(fileName,suff);end %By defult .mat file assumed

%%% Set default flag
flag=0;

%%% Check if file exists
if exist(fileName,'file')==2
    flag=flag+1;
    
    if ~isempty(varName) && strcmp( suff,'.mat')
        variableInfo = who('-file',fileName);
        flag = flag + single(ismember(varName, variableInfo));
    end
    
end


