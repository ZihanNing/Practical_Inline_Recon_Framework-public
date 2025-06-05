function [user,versCode,versSave,pathRF,pathSt,pathPr,pathNt]=versRecCode(rootOu)

%VERSRECCODE   Returns the current version, different paths and options for
%the reconstruction
%   [USER,VERSCODE,VERSSAVE,PATHRF,PATHST]=VERSRECCODE
%   * USER is the user executing the code
%   * VERS is the version of the code
%   * VERSSAVE is the matlab 
%   * PATHRF is the absolute path to the ReconFrame version
%   * PATHST is the path to the specific studies information
%   * PATHPR is the path to project information files
%   * PATHNT is the path to pretrained neural networks
%

%MODIFY THIS WITH MOST CONVENIENT TO YOUR NEEDS (LINUX USERS, IN WINDOWS FULL PATHS ARE USED)
user='lcg13';
versCode='Release07';
versSave='-v7.3';%Faster than -v7.3 due to lack of compression but limited to 2GB files
pathSt=fullfile(fileparts(mfilename('fullpath')),'..','Data');
%pathSt=fullfile(filesep,'home',user,'Work',strcat('DefinitiveImplementation',versCode),'Data');
%If using MRecon-3.0.460 we get 5x quicker reading and 2x less storage in dumped headers but without the GoalC information

%MODIFY THIS WITH YOUR MRECON/STUDIES PATH
%pathRF=fullfile(filesep,'home',user,'Work','MRecon-3.0.460');
pathRFV{1}=fullfile(fileparts(mfilename('fullpath')),'..','..','MRecon-3.0.515');
pathRFV{2}=fullfile(filesep,'home',user,'Work','MRecon-3.0.515');
pathRFV{3}=fullfile(filesep,'usr','local','MATLAB','MRecon-3.0.515');
%pathRF=fullfile(filesep,'home',user,'Work','MRecon-3.0.537');

pathRF=[];
for s=1:length(pathRFV)
    if exist(pathRFV{s},'dir')
        pathRF=pathRFV{s};
        break;
    end
end

%WE MAKE THIS SEARCH INDEPENDENT OF THE OUTPUT FOLDER
%if nargin<1 || isempty(rootOu);pathPrOrV{1}=fullfile(filesep,'home',user,'Data','pnrawDe');
%else pathPrOrV{1}=fullfile(rootOu,'..');
%end
pathPrOrV{1}=fullfile(filesep,'pnraw01','dhcp-reconstructions');
pathPrOrV{2}=fullfile(filesep,'home','lcg13','Data','pnrawDe');
pathPrOr=[];
for s=1:length(pathPrOrV)
    if exist(pathPrOrV{s},'dir') && length(dir(pathPrOrV{s}))>2%More than just . and ..
        pathPrOr=pathPrOrV{s};
        break;
    end
end

%assert(~isempty(pathRF),'Unable to find ReconFrame folder');

pathPr{1}=fullfile(pathPrOr,'StudiesReplicated');%Neonatal
pathPr{2}=fullfile(pathPrOr,'StudiesReplicatee');%Fetal
pathPr{3}=fullfile(pathPrOr,'Studies','DanaRelease06.txt');%Toddlers
pathPr{4}=fullfile(pathPrOr,'Studies','BrainInjuryRelease03.txt');%Brain injury
pathPr{5}=fullfile(pathPrOr,'Studies','EpilepsyRelease07.txt');%Toddlers

pathNt=fullfile(filesep,'home',user,'Work','TrainedNetworks');

%pathRF=fullfile('c:\Users\310187552\Documents\MATLAB\MRRecon');%Torben

%pathRF=fullfile(filesep,'usr','local','MATLAB','MRecon-3.0.515');

%See below the results of an experiment on beastie01 (repeated later for perinatal130-pc, see relative times) comparing
%pathRF=sprintf('/home/%s/%s/MRecon-3.0.460',user,matlabRoot); with
%pathRF='/usr/local/MATLAB/MRecon-3.0.515'; for some phantom
%(neonatal/fetal geometry/fetal dual echo) data.
%
%
% %pathRF=sprintf('/home/%s/%s/MRecon-3.0.460',user,matlabRoot);
% 
% Parsing raw data folder /home/lcg13/Data/rawSource/archive-rawdata/archive-nnu/2015_10_20/dH_64800
% The destination folder already existed, overwritting parsed content
% Size of .mat files: 1.17 Gb. Ratio of compression over .lab: 14.23 %
% Time parsing raw data folder: 137.953 s - 94.748 s
% Parsing raw data folder /home/lcg13/Data/rawSource/archive-rawdata/archive-nnu/2017_01_18/ge_193200
% The destination folder already existed, overwritting parsed content
% Size of .mat files: 512.10 Mb. Ratio of compression over .lab: 26.03 %
% Time parsing raw data folder: 58.901 s - 49.403 s
% Parsing raw data folder /home/lcg13/Data/rawSource/pnraw/raw-nnu/2017_08_29/ba_500
% The destination folder already existed, overwritting parsed content
% No scan technique stored for series ba_29082017_1637228_5_2_pud2b0mapshimmo1V4
% Size of .mat files: 166.32 Mb. Ratio of compression over .lab: 10.76 %
% Time parsing raw data folder: 33.540 s - 28.028 s
% 
% 
% %pathRF='/usr/local/MATLAB/MRecon-3.0.515';
% 
% Parsing raw data folder /home/lcg13/Data/rawSource/archive-rawdata/archive-nnu/2015_10_20/dH_64800
% The destination folder already existed, overwritting parsed content
% Size of .mat files: 1.17 Gb. Ratio of compression over .lab: 14.23 %
% Time parsing raw data folder: 136.887 s - 103.984 s
% Parsing raw data folder /home/lcg13/Data/rawSource/archive-rawdata/archive-nnu/2017_01_18/ge_193200
% The destination folder already existed, overwritting parsed content
% Size of .mat files: 594.78 Mb. Ratio of compression over .lab: 41.75 %
% Time parsing raw data folder: 245.280 s - 231.076 s
% Parsing raw data folder /home/lcg13/Data/rawSource/pnraw/raw-nnu/2017_08_29/ba_500
% The destination folder already existed, overwritting parsed content
% Size of .mat files: 255.08 Mb. Ratio of compression over .lab: 16.42 %
% Time parsing raw data folder: 112.358 s - 89.416 s
