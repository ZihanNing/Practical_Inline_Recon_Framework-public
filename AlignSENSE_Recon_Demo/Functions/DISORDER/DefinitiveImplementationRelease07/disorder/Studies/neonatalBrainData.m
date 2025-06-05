%path{1}='2015_06_19/PO_36300';%107-Processed for Judit to see SNR helmet effect in fMRI
%path{2}='2016_07_29/CH_144900';%498-Processed for Judit to see SNR helmet effect in fMRI

%path{1}='2015_03_24/LU_18600';%Neonatal subject with intermediate quality as reported in Matteo's paper

%path{1}='2018_05_14/MA_11130';%Processed for Tanya, not acquired with dHCP patch

%path{3}='2018_09_03/HE_43130';%Neonatal case with spikes, to study the MPRAGE
%path{1}='2018_08_21/AB_39730';%Neonatal case with spikes, to study
%path{2}='2018_08_24/WA_41130';%Neonatal case with spikes, to study

%path{1}='2018_09_06/FE_44030';
%path{1}='2018_09_06/HE_43930';

%path{1}='2017_01_09/MA_190801';%Neonatal case with corrupted information in series 5 

%Old diffusion tests to set up the protocol
%path{1}='2014_06_12/TU_1601';%b-val:{[9 10 11 12 13]}
%path{2}='20140325/BH';%b-val:{[12 13 14 15 16]}
%path{3}='2014_09_10/SU_18303';%b-val:{[6 7 8 9 10]}
%path{4}='2014_09_30/AL_22402';%b-val:{[6 7 8 10 11]}

%path{1}='2019_01_03/HA_70730';%Brain injury Joubert's disease

%path{1}='2018_12_04/PE_66031';%Issues with the coils probably---PROBLEMS
%path{2}='2018_12_17/LE_68830';%Issues with the coils probably---PROBLEMS
%path{3}='2018_12_21/WI_70230';%Perhaps issues with the coils---NOT OBSERVED PROBLEMS BUT PROBABLY THERE ARE SOME, GOOD TO RECHECK
%path{4}='2019_01_10/AD_72130';%Perhaps issues with the coils---PROBLEMS
%path{5}='2019_01_10/AD_72230';
%path{6}='2019_01_10/AD_72330';
%path{7}='2019_01_11/TA_72531';
%path{8}='2019_01_14/DE_72930';

%path{1}='2019_01_25/BA_76130';%Dramatic failure with new coil (that from December)

%path{1}='2019_01_18/co_74430';%New test with the returned coil
%path{2}='2019_01_18/ST_74431';%New neonatal scan

%path{1}='2019_01_17/BO_74030';%Potential issue with fMRI data- coil detuning observed

%path{1}='2019_01_19/HA_74630';%8-year old
%path{2}='2019_01_19/MA_74530';%8-year old
%path{1}='2019_01_22/LE_75031';%BIG ISSUE WITH RETURNED COIL


%path{5}='2019_01_11/RU_72430';%Check in a fetal case-nothing

%path{1}='2019_02_05/KW_78730';%Case where detuning has first appeared in another channel, 19, instead of 29
%path{2}='2019_02_08/BA_79430';%Case where detuning has first appeared in channel 24.

%ZEBRA STUDY JANA:
%path{1}='2019_03_07/CH_86031';%Scan number 18 is the zebra data
%path{1}='2019_03_13/MA_87730';%Scans 7-10 were to be corrected

%NEONATAL DISTORTION CORRECTION
%path{1}='2015_11_26/SU_78700';%Interesting fMRI case to check motion and distortion correction in neonates-%CC00236XX14 78700
%CASES SPOTTED BY JUDIT:
%path{1}='2015_07_30/LA_45100';%Premature
%path{2}='2017_05_23/MI_222000';%Hematoma
%path{3}='2015_08_10/BA_47600';%Deformed skull
%path{4}='2015_10_01/RO_59500';%Big motion
%path{5}='2016_03_01/ST_104100';%Big CSF
%path{6}='2015_06_30/UG_38200';%Normal
%Note this is the way to run fix: http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/centos6/fsl_507/doc/wiki/FIX.html
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/neofmri_2nd_release_rerun2/sub-CC00157XX09/ses-51900/fix/filtered_func_data.nii.gz-Preprocessed data
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/neofmri_2nd_release_rerun2/sub-CC00157XX09/ses-51900/fix/filtered_func_data_clean.nii.gz-Fix-cleaned data
%Hints on the pipeline on:
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/seanf/dhcp-pipeline/README.md
%Graphically:
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/seanf/dhcp-pipeline/dhcp/schematics/pipeline_dhcp.dot.png
%The pipeline in a script:
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/seanf/scripts/run.sh
%Actually motion correction is in script:
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/neofmri_2nd_release_rerun2/sub-CC00236XX14/ses-78700/logs/mcdc.log
%Motion and distortion correction results in:
%/home/lcg13/Data/pnrawOrRelease03/dhcp-pipeline-data/oxford/dhcp-results/neofmri_2nd_release_rerun2/sub-CC00236XX14/ses-78700/mcdc/func_mcdc.nii.gz-Check whether this matches filtered_func_data.nii.gz

%MORE SPIKES:
%path{1}='2019_05_13/CO_99330';%VERY STRONG IN DIFFUSION BUT APPEAR
%UNSTRUCTURED-10 is the B1

%COIL 2-CARDIAC STUDY-CHECKING FOR DETUNING IN MPRAGE?
%path{1}='2019_06_04/VU_104331';

%path{1}='2015_07_01/CH_38800';%TO CHECK PARSING OF SHIM_B0 AT THE END

%path{1}='2015_05_07/CO_27000';%dHCP case to test reconstruction%CC00071XX06
%path{1}='2015_05_12/SH_28000';%dHCP case to test reconstruction challenging fMRI-better now!%CC00074XX09
%path{1}='2015_05_13/GE_28400';%dHCP case to test artifact in diffusion%CC00075XX10

%path{1}='2018_08_13/PE_37430';%sub-CC00928XX21 ses-37430
%path{2}='2015_07_20/KA_43100';%sub-CC00126XX11 ses-43100
%path{3}='2016_07_29/CH_144900';%sub-CC00498XX21 ses-144900
%path{4}='2016_10_21/CR_170600';%sub-CC00569XX17 ses-170600
%path{5}='2016_12_20/TA_188800';%sub-CC00592XX16 ses-188800
%path{6}='2015_05_13/GE_28400';%sub-CC00075XX10 ses-28400

%%FOR RITA, MOTION FREE
%%path{1}='2015_03_03/AL_13300';%It has another scan at the end
%path{1}='2015_05_18/PI_29500';
%path{2}='2015_05_21/DI_30500';
%path{3}='2015_05_28/CH_31800';
%path{4}='2015_05_29/MI_31801';
%path{5}='2015_06_01/CH_32100';

% path{1}='2015_06_10/DI_33702';
% path{1}='2015_06_10/MA_33700';
% path{1}='2015_06_18/TA_36100';
% path{1}='2015_06_24/GI_37001';
% path{1}='2015_07_02/SA_39400';
% path{1}='2015_07_30/LA_45000';
% path{1}='2015_08_13/AK_48500';
% path{1}='2015_08_24/VA_50700';
% path{1}='2015_10_08/LA_61000';
% path{1}='2015_11_27/GB_79300';
% path{1}='2015_12_10/KI_83000';
% path{1}='2016_01_14/LA_89400';
% path{1}='2016_02_11/CA_97401';
% path{1}='2016_03_31/OZ_113001';
% path{1}='2016_04_11/WI_115700';

%%MOTION CORRUPTED
%path{6}='2015_06_15/NA_35000';
%path{7}='2015_06_17/DO_35801';
%path{8}='2015_12_21/TY_85100';
%path{9}='2016_03_17/JI_109200';
%path{10}='2016_04_07/DO_114500';

%path{1}='2019_10_02/KL_130330';

%path{1}='2019_10_16/NO_132530';%Neonate acquired with the adult coil

%DISORDER FIRST TEST ON NEONATES
%path{1}='2017_02_01/PE_197500';
%path{2}='2017_02_15/VA_201400';

%%DISORDER QMRI-NEONATES
%path{1}='2019_11_18/BO_139230';
%path{2}='2019_11_18/HO_139231';

%%DISORDER QMRI-NEONATES
%path{1}='2019_11_27/CA_140730';

%path{1}='2019_11_28/VE_141130';%Neonate coil issues, looking good

%MOTION CORRECTION
%path{1}='2015_02_27/NG_12501';%Neonatal case used for T2MS corrections in the paper%CC00060XX03
%path{2}='2015_09_09/CO_54400';%Neonatal case used for T1MS corrections in the paper%CC00135AN12

%DAAN RECONSTRUCTIONS
%path{1}='2016_07_29/CH_144900';

%FMRI TESTS
%path{1}='2015_03_03/AL_13300'; %CC00061XX04---HERE THERE IS A PROBLEM OF INHOMOGENEITY, THIS IS CAUSED BY CHANNEL 18
path{2}='2015_05_07/CO_27000'; %CC00071XX06---CASE 23
path{3}='2015_05_12/SH_28000'; %CC00074XX09
path{4}='2015_06_02/TO_32300'; %CC00091XX10
path{5}='2015_06_17/AH_35800'; %CC00104XX05
path{6}='2015_07_21/RA_43500'; %CC00129AN14
path{7}='2015_11_18/GI_76601'; %CC00227XX13
%path{8}='2016_03_17/JI_109200';%CC00346XX17
