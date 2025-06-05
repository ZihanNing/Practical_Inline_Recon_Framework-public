function rec=espiritAlgorithm(rec)

%ESPIRITALGORITHM   Sets the default parameters for ESPIRiT
%   REC=ESPIRITALGORITHM(REC)
%   * REC is a reconstruction structure without the algorithmic parameters. 
%   At this stage it may contain the naming information (.Names), the 
%   status of the reconstruction (.Fail), the .lab information (.Par), the 
%   fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   ** REC is a reconstruction structure with filled algorithmic 
%   information (.Alg)
%


%ALGORITHM
if ~isfield(rec,'Alg');rec.Alg=[];end

%DATA
if ~isfield(rec,'y');error('Data not present. Please provide data in field rec.y');end
if ~isfield(rec,'S');error('Sensitivities not present. Please provide sensitivities in field rec.S');end
if ~isfield(rec,'M');error('Mask not present. Please provide mask in the field rec.M');end

%OUTPUT FILE
if ~isfield(rec,'Names');rec.Names=[];end
if ~isfield(rec.Names,'pathOu') || isempty(rec.Names.pathOu);error('No output path provided. Please provide one in field rec.Names.pathOu');end
if ~isfield(rec.Names,'Name') || isempty(rec.Names.Name);error('No output file provided. Please provide one in field rec.Names.Name');end

%TRAJECTORIES
if ~isfield(rec,'Assign');rec.Assign=[];end
if ~isfield(rec.Assign,'z') || ~iscell(rec.Assign.z) || length(rec.Assign.z)<3 || isempty(rec.Assign.z{2}) || isempty(rec.Assign.z{3}) || numel(rec.Assign.z{2})~=numel(rec.Assign.z{3});error('Trajectories not set. These should go into a cell rec.Assign.z with second element for first PE (row vector) and third element for second PE (page vector). Both elements should have the same number of entries, which should correspond to the number of profiles');end

%SCAN INFO
if ~isfield(rec,'Par');rec.Par=[];end
if ~isfield(rec.Par,'Labels');rec.Par.Labels=[];end
if ~isfield(rec.Par.Labels,'TFEfactor')
    fprintf('Number of samples per effective shot (rec.Par.Labels.TFEfactor) not set. Default is used (only valid for steady-state)\n');
    rec.Par.Labels.TFEfactor=length(rec.Assign.z{2});
end    
if ~isfield(rec.Par.Labels,'ZReconLength')
    fprintf('Number of effective shots (rec.Par.Labels.ZReconLength) not set. Default is used (only valid for steady-state)\n');
    rec.Par.Labels.ZReconLength=1;
end

%RESOLUTION
if ~isfield(rec,'Enc');rec.Enc=[];end
if ~isfield(rec.Enc,'AcqVoxelSize')
    fprintf('Resolution (rec.Enc.AcqVoxelSize) not set. Default 1mm isotropic is used\n');
    rec.Enc.AcqVoxelSize=ones(1,3);%Resolution to 1mm
end
NY=size(rec.y);
if ~isfield(rec.Enc,'kRange');rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};end %YB: factor 2 comes from the oversampling in the readout direction (1st dimension here)
if ~isfield(rec.Enc,'FOVSize');rec.Enc.FOVSize=NY(1:3);end
if ~isfield(rec.Enc,'AcqSize');rec.Enc.AcqSize=[2*NY(1) NY(2:3)];end

%GEOMETRY
if ~isfield(rec.Par,'Mine');rec.Par.Mine=[];end
if ~isfield(rec.Par.Mine,'Asca');rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);end
if ~isfield(rec.Par.Mine,'Arot');rec.Par.Mine.Arot=eye(4);end
if ~isfield(rec.Par.Mine,'Atra');rec.Par.Mine.Atra=eye(4);end
if ~isfield(rec.Par.Mine,'MTT');rec.Par.Mine.MTT=inv([-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]);end
if ~isfield(rec.Par.Mine,'APhiRec');rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;end

%SCAN GEOMETRY
% if ~isfield(rec.Par,'Scan');rec.Par.Scan=[];end
% if ~isfield(rec.Par.Scan,'MPS');rec.Par.Scan.MPS='HF LR PA';end
% if ~isfield(rec.Par.Labels,'FatShiftDir');rec.Par.Labels.FatShiftDir='H';end
% if ~isfield(rec.Par.Labels,'FoldOverDir');rec.Par.Labels.FoldOverDir='PA';end%This is the quick/first phase encode

%SCAN PARAMETERS
if ~isfield(rec.Par.Labels,'ScanDuration')
    fprintf('Scan duration not set. Default 1000s is used, temporal plots may not look well\n');
    rec.Par.Labels.ScanDuration=1000;
end
if ~isfield(rec.Par.Labels,'RepetitionTime')
    fprintf('Repetition time not set. Default 1s is used, temporal plots may not look well\n');    
    rec.Par.Labels.RepetitionTime=1;%1 second by default
end
if ~isfield(rec.Par.Labels,'TE')
    fprintf('Echo time not set. Default 1ms is used, temporal plots may not look well\n');    
    rec.Par.Labels.TE=0.001;%1 milisecond by default
end
if ~isfield(rec.Par.Labels,'FlipAngle');rec.Par.Labels.FlipAngle=0;end
if ~isfield(rec.Par.Scan,'Technique');rec.Par.Scan.Technique='T1TFE';end
if ~isfield(rec.Par.Mine,'Modal');rec.Par.Mine.Modal=7;end%Reconstruction modality - 7 refers to 

%ALGORITHM (SPECIFIC)

%DATA TYPE INFORMATION
if ~isfield(rec,'Plan');rec.Plan=[];end
if ~isfield(rec.Plan,'Dims');rec.Plan.Dims={'size','ky','kz','chan','dyn','card','echo','loca','mix','extr1','extr2','aver'};rec.Plan.NDims=length(rec.Plan.Dims);end
%1-> Raw / 2-> Empty / 3-> EPIRead / 4-> Navigator / 5-> Noise / 6-> Spectra
%7-> Sensitivities / 8-> Masks / 9-> G-facts / 10-> Chi2-maps / 11-> B0 / 12-> Reconstruction
%13-> Undistorted / 14-> Per-volume alignment / 15-> Per-excitation alignment / 16-> Volumetric alignment / 17-> Motion transform / 18-> Filtered reconstruction
%19-> Noise level / 20-> Number of components / 21-> Asymptotic error / 22-> Residuals / 23-> Frequency stabilization data / 24-> Image after frequency stabilization
%25-> Tracking information / 26->Tracked data / 27 -> SensitivityEigenMaps / 28 -> Chi2-maps in material coordinates /  29 -> Linear B0 varying with rotation
if ~isfield(rec.Plan,'Types');rec.Plan.Types={'z','','P','C','N','y', ...
            'S','M','G','X','B','x', ...
            'u','v','e','d','T','r', ...
            's','p','a','n','F','w', ...
            't','b','W','E','D'};end
if ~isfield(rec.Plan,'TypeNames');rec.Plan.TypeNames={'Ra','','Ny','Na','Ga','Sp', ...
                'Se','Ma','No','Ch','B0','Aq', ...
                'Un','Vo','Ex','Di','Tr','Re', ...
                'Si','Co','Ae','Er','Fr','Ws', ...
                'Fo','Vt','Ei','Ec','Dc'};rec.Plan.NTypes=length(rec.Plan.TypeNames);end
if ~isfield(rec.Plan,'CorrTypes');rec.Plan.CorrTypes={'random_phase','pda_fac','pda_index','meas_phase','sign'};rec.Plan.NCorrs=length(rec.Plan.CorrTypes);end

%ACTUALLY USED INFORMATION
if ~isfield(rec,'Dyn');rec.Dyn=[];end

%ERROR HANDLING
if ~isfield(rec,'Fail');rec.Fail=0;end%Flag to abort

%GPU
if ~isfield(rec.Dyn,'BlockGPU');rec.Dyn.BlockGPU=1-double(gpuDeviceCount>0);end%Flag to block GPU computations
if ~isfield(rec.Dyn,'MaxMem');rec.Dyn.MaxMem=[6e6 2e6 1e6];end%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

%OTHER
if ~isfield(rec.Plan,'SuffOu');rec.Plan.SuffOu='';end%Suffix to add to output data
if ~isfield(rec.Plan,'Suff');rec.Plan.Suff='';end%Suffix to add to images
if ~isfield(rec.Par.Mine,'Proce');rec.Par.Mine.Proce=0;end%Flag to indicate that this is a reconstruction problem instead of a preprocessing problem
if ~isfield(rec.Par.Mine,'Signs');rec.Par.Mine.Signs=[];end%For rotating PEs, not used
if ~isfield(rec.Par.Mine,'Nat');rec.Par.Mine.Nat=1;end%Native PE in case of rotating PEs, not used
if ~isfield(rec.Par,'Encoding');rec.Par.Encoding=[];end%Additional encoding information, not used
