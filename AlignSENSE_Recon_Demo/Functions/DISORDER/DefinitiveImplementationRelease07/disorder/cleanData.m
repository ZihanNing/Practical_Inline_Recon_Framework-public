function cleanData(path,modal,rootOu)

%addpath(genpath(fileparts(mfilename('fullpath'))));


%path='2017_08_29/ba_500';
%path='2017_01_18/ge_193200';
%path='2015_10_20/dH_64800';
%path='2017_05_23/ph_222300';

%CLEANDATA   Cleans existing reconstructions
%   PARSELABFOLDER(PATH,{MODAL},{ROOTOU})
%   * PATH is the relative path to raw and parsed data, for instance 
%   2014_04_08/LO_10203
%   * {MODAL} is a cell array / a vector / empty (default) for the set of
%   modalities of interest. Empty corresponds to all contemplated
%   modalities
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions)
%

%SET DEFAULT PARAMETERS AND DETECT THE RAW AND PARSED DATA FOLDERS
[user,versCode]=versRecCode;

if ~exist('modal','var');modal=[];end
if ~exist('rootOu','var');rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end
pathOu=fullfile(rootOu,path);
if isempty(modal) || ~isnumeric(modal);modal=modal2Numbe(modal);end

for n=modal
    pathModa=fullfile(pathOu,numbe2Modal(n));
    if exist(pathModa,'dir');rmdir(pathModa,'s');end
end
