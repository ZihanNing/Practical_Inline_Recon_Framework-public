

function [] = transferData(strContain, type, sourceMach, destMach, sourceDataDir,destDataDir,userName)

%TRANFERDATA  transfers data from one machine to the other, with the
%option to specify directories and types of data to be transfered. Note
%that this uses the Linux command line so might not work on Windows and
%might require manual password entry (more than once if multiple types
%copied). For the latter issue, consider setting up password keys on the
%between which to copy files. 
%   []=TRANFERDATA(STRCONTAIN, {TYPE},{SOURCEMACH},{DESTMACH},{SOURCEDATADIR},{DESTDATADIR},{USERNAME} )
%   * STRDIR is a string to is contained in the directory (relative to the dataDir) to be transfered.
%   * {TYPE} is a cell structure with strings of data types to be transfered. Defaults to all.
%   * {SOURCEMACH} is the machine were to copy data from.
%   * {DESTMACH} is the machine were to copy data to.
%   * {SOURCEDATADIR} is the directory where the data is stored on the source machine.
%   * {DESTDATADIR} is the directory where the data is stored on the destiny machine.
%   * {USERNAME} is the username under which login to transfer the dat.
%

allTypes = {'Logging','Re-Se','Re-F0','Re-B1','An-Ve','Re-Se_Sn','An-Ve_Sn','ISD-PACS','rawmat','rawdat','rawmeta'};

if nargin<2 || isempty(type); type = allTypes;end%By default,copy everything
if nargin<3 || isempty(sourceMach); sourceMach = '';end%'gpubeastie05'
if nargin<4 || isempty(destMach); destMach = 'perinatal005-pc';end
if nargin<5 || isempty(sourceDataDir); sourceDataDir = '/home/ybr19/Data/';end
if nargin<6 || isempty(destDataDir); destDataDir = sourceDataDir;end
if nargin<7 || isempty(userName); userName = 'ybr19';end

addColonSource = ~strcmp(sourceMach,'');
addColonDest = ~strcmp(destMach,'');

if addColonSource; colonSource = ':';nameSource=strcat(userName,'@');else; colonSource='';nameSource='';end
if addColonDest; colonDest = ':';nameDest=strcat(userName,'@');else; colonDest='';nameDest='';end

flag = cellfun( @(x) strcmp(x,'raw'), type); 
if any(flag); type(flag==1)=[]; type(end+1:end+2) = {'rawmat','rawdat'};end

%%% Extract the directory to transfer
dirInfo = dir(fullfile( sourceDataDir) );
dirInfo( ~[dirInfo.isdir]) = []; %Keep directories

Names = {dirInfo.name};

flag = contains( Names, strContain) ;
if multDimSum(flag)>1
    warning('transferData:: Multiple directories detected. Aborted data transfer.');
    return;
elseif multDimSum(flag)<1
    warning('transferData:: No directories detected. Aborted data transfer.');
    return;
end

idx = find(flag);
dirName = Names{idx(1)}(1:end);%Only first to transfer

%%% Run over types and transfer
for i = 1:length(type)
    fprintf('\nTransferring %s:\n',type{i});
    if strcmp(type{i},'rawmat')%Matlab converted raw data
        command = sprintf('scp -r %s%s%s%s%s/*.mat %s%s%s%s%s/',...
                           nameSource,sourceMach,colonSource,sourceDataDir,dirName,...
                           nameDest,destMach,colonDest,destDataDir,dirName);
        dos(command);
        
    elseif strcmp(type{i},'rawdat') %Siemens raw data
        command = sprintf('scp -r %s%s%s%s%s/*.dat %s%s%s%s%s/',...
                           nameSource,sourceMach,colonSource,sourceDataDir,dirName,...
                           nameDest,destMach,colonDest,destDataDir,dirName);
        dos(command);
        
    elseif strcmp(type{i},'rawmeta')%JSON metedata
        command = sprintf('scp -r %s%s%s%s%s/*.json %s%s%s%s%s/',...
                           nameSource,sourceMach,colonSource,sourceDataDir,dirName,...
                           nameDest,destMach,colonDest,destDataDir,dirName);
        dos(command);
        
    else
        command = sprintf('scp -r %s%s%s%s%s/%s %s%s%s%s%s/',...
                           nameSource,sourceMach,colonSource,sourceDataDir,dirName,type{i},...
                           nameDest,destMach,colonDest,destDataDir,dirName);%No type{i} here as need to prive parent dir
        dos(command);
    end
    
end

