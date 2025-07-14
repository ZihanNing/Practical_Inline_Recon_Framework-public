[~,hostname]=system('hostname');hostname=strtrim(lower(hostname));

% raw data path
pathData = '/home/gadgetron/Gadgetron_Parallel_Framework';
caseIn='raw';
bodyIn='';
dataIn{1}='RAW_SWI_forSENSE_integratedACS.dat';
surfIn=dataIn; % use embedded ACS line for reconstruction; could also be replaced by a reference dataset for coil sensitivity map estimation

if ~iscell(surfIn);surfIn=repmat({surfIn},1,length(dataIn));end
caseIn=strcat(pathData,filesep,caseIn);