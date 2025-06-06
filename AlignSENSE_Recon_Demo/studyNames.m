
pathIn{1}='/home/gadgetron/Gadgetron_Parallel_Framework/raw';%Where on the server the data is (temporally stored) NEEDS TO HAVE THE / AT THE END

%% ACQUISITIONS
fileIn{1}{1}='meas_MID00107_FID16038_T1_MPRAGE_DISORDER_twinsACS'; 
% fileIn{1}{2}='meas_MID00054_FID02071_SWI_DISORDER_fullsampled_10_10_BW500_pose2'; 
% fileIn{1}{3}='meas_MID00055_FID02072_SWI_DISORDER_fullsampled_10_10_BW500_pose3'; 

%% REFERENCE SCANS
refIn{1}{1} = '';
refBIn{1}{1} = '';
% refIn{1}{1} = 'meas_MID00029_FID64029_MPRAGE_DISORDER_p2_RandSeed22LIP_REF';

%% DATA CONVERSION SPECIFICATIONS - normally no need to touch this (only in case you accidentally acquired data at a too high resolution or FOV)
RDesired = fillCell(fileIn,[]); %e.g.[ 1.2 1.3]
