[~,hostname]=system('hostname');hostname=strtrim(lower(hostname));

% ZIHAN
pathData = '/home/zn23/Gadgetron_Parallel_Framework';
caseIn='raw';
bodyIn='';
dataIn{1}='meas_MID00031_FID85168_SWI_UKBB_Gadget_20250522_002252.dat';
surfIn=dataIn;

if ~iscell(surfIn);surfIn=repmat({surfIn},1,length(dataIn));end
caseIn=strcat(pathData,filesep,caseIn);

%fMRI_2023_07_04