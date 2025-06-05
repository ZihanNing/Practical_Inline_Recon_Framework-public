
function rec = ISMRMRD2Rec(fileName, supportReadout, writeNIIFlag, writeRAW, pathOu, isPilotTone)

%ISMRM2REC Converts ISMRM raw data format (ISMRM-RD) to a reconstruction structure.
%   * FILENAME is the filename of the raw data file.
%   * {SUPPORTREADOUT} is the support in the readout direction for high resolution scans.
%   * {WRITENIIFLAF} a flag whether to save the raw data into a .mat structure.
%   * {SHIMDIR} is directory with additional shim information of the scanner
%   ** REC is a reconstruction structure
%

if nargin<2 || isempty(supportReadout); supportReadout=[];end
if nargin<3 || isempty(writeNIIFlag); writeNIIFlag=[];end
if nargin<4 || isempty(writeRAW); writeRAW=[];end
if nargin<5 || isempty(pathOu); pathOu=[];end
if nargin<6 || isempty(isPilotTone); isPilotTone=0;end

gpu=(gpuDeviceCount>0 && ~blockGPU);

%% OPEN H5 ISMRMRD FILE (need master from https://github.com/ismrmrd/ismrmrd)
if exist(strcat(fileName,'.h5'), 'file')
    h5File = ismrmrd.Dataset(strcat(fileName,'.h5'));
else
    error(['File ' strcat(fileName,'.h5') ' does not exist. Please generate it using the ISMRMRD converter.'])
end
h5disp(strcat(fileName,'.h5'))
TEMP = h5read(strcat(fileName,'.h5'), '/dataset/data');

%% READ DATA
D = h5File.readAcquisition();%Acquisition
hdr = ismrmrd.xml.deserialize(h5File.readxml);%Header

%%% ASSIGN NOISE
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan>1; noise = D.select(1:firstScan-1); else; noise = []; end

%%% ASSIGN MEASUREMENTS
meas = D.select(firstScan:D.getNumber); 
clear D

% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;
%NEnc = 
%NRec = 

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
%FOVEnc
% rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
% rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
% rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

if ~isempty(supportReadout);vr=round(rec_Nx*supportReadout(1)):round(rec_Nx*supportReadout(2));else; vr=[];end

if isfield(hdr.encoding.encodingLimits.slice,'maximum'); nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;else;  nSlices = 1;end
if isfield(hdr.acquisitionSystemInformation,'receiverChannels'); nCoils = hdr.acquisitionSystemInformation.receiverChannels;else;  nCoils = 1;end
if isfield(hdr.encoding.encodingLimits.repetition,'maximum'); nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;else;  nReps = 1;end
if isfield(hdr.encoding.encodingLimits.contrast,'maximum'); nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 +1;else;  nContrasts = 1;end

%% EXTRACT COIL K-SPACE DATA
%%% NOISE
rec.N=[];
for p = 1:length(noise.data)
    rec.N = cat(3,rec.N,noise.data{p});
end
clear noise

rec.N = permute(rec.N, [1 3 4 2]);
if isfield(rec,'N') && gpu;rec.N=gpuArray(rec.N);end
for l=1
    rec.N= fftshiftGPU(rec.N,l); 
    rec.N= fftGPU(rec.N,l);     %*NY(l)  YB: this has once given inf as number close to realmax('single') 
    rec.N= ifftshiftGPU(rec.N,l);
end
Nnoise = size(rec.N);Nnoise(1) = rec_Nx;
if isfield(rec,'N');rec.N=resampling(rec.N,Nnoise,2);end

%%% ACQUISITION
y = [];yTemp1=[];yTemp2=[];
for rep = 1:nReps
    for contrast = 1:nContrasts
        for slice = 1:nSlices
            % Initialize the K-space storage array
            K = zeros(enc_Nx, enc_Ny, enc_Nz, nCoils);
            
            % Select the appropriate measurements from the data
            acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
                        & (meas.head.idx.repetition==(rep-1)) ...
                        & (meas.head.idx.slice==(slice-1)));
            for p = 1:length(acqs)
                ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                K(:,ky,kz,:) = permute( meas.data{acqs(p)},[1 3 4 2]);
            end
            if gpu;K=gpuArray(K);end
            K = single(K);
            
            % To image domain in RO
            for l=1
                K= fftshiftGPU(K,l); 
                K = fftGPU(K,l);  %YB: was ifft originally
                K = ifftshiftGPU(K,l);
            end
            
            % Remove over-sampling
            NK= size(K);NK(1) = rec_Nx;
            K=resampling(K,NK,2);
            
            % To image domain in PE1 and PE2
            for l=2:3
                K= fftshiftGPU(K,l); 
                K = fftGPU(K,l);  %YB: was ifft originally
                K = ifftshiftGPU(K,l);
            end
            
            yTemp1 = cat(6,yTemp1,K); 
        end
        yTemp2 = cat(7,yTemp2,yTemp1);
        yTemp1=[];
    end
    y = cat(5,y,yTemp2);
    yTemp2=[];
end

rec.y = y; y =[];
if isfield(rec,'N');rec.y = standardizeCoils(rec.y,rec.N);end
if gpu
    rec.y = gather(rec.y);
    if isfield(rec,'N'); rec.N = gather(rec.N);end
end

%% GEOMETRY COMPUTATION
%SCALING
rec.Enc.AcqVoxelSize = single( [ enc_FOVx/enc_Nx enc_FOVy/enc_Ny enc_FOVz/enc_Nz]);
rec.Par.Mine.Asca = diag([ rec.Enc.AcqVoxelSize 1]);

%ROTATION
rec.Par.Mine.Arot = single( eye(4));
rec.Par.Mine.Arot(1:3,1:3)=single( cat(2, meas.head.read_dir(:,1), meas.head.phase_dir(:,1), meas.head.slice_dir(:,1))); %RO-PE-SL to PCS

%TRANSLATION
rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = single(meas.head.position(1:3,1));%RO-PE-SL to PCS

% account for fact that  meas.head.position is referred to the center of FOV, not the first element in the array
NY = size(rec.y); NY=NY(1:3);
orig = ((NY+0)/2 - 0 - [0 0 .5])'; orig=orig(1:3);
if ~isempty(vr); orig(2)=orig(2)-(vr(1)-1); end%YB:orig(2) since this is readout/only vr(1)-1 elements removed from array/minus sign since should be addition and below already negative sign

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)- rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;
                    
%COMBINED MATRIX
rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

%MOVE TO RAS
rec.Par.Mine.PCS2RAS = diag([-1 -1 1 1]);%For HeadFirstSupine. Can be generalised from hdr.measurementInformation.patientPosition
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From RO-PE-SL to RAS
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;

%DEDUCE ACQUISITION ORDER
[~,MT_acqOrder] = mapNIIGeom([], rec.Par.Mine.APhiRec, 'permute', [2 1 3], NY);%acquisitionOrder takes MT as input that is defined in the PRS to RAS
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane ] = acquisitionOrder(MT_acqOrder);clear MT_acqOrder;
rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
fprintf('\nSlice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
fprintf('Readout: %s\n', rec.Par.Scan.Mine.RO );
fprintf('1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
fprintf('2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO

%% SEQUENCE INFORMATION
rec.Par.Labels.TFEfactor=length(meas.head.idx.kspace_encode_step_1);%This should be changed for shot-based sequences       
rec.Par.Labels.ZReconLength=1;

rec.Par.Labels.RepetitionTime=hdr.sequenceParameters.TR;%In ms
rec.Par.Labels.TE=hdr.sequenceParameters.TE;%In ms
rec.Par.Labels.FlipAngle=hdr.sequenceParameters.flipAngle_deg;%In degrees
rec.Par.Labels.ScanDuration= rec.Par.Labels.RepetitionTime * length(meas.head.idx.kspace_encode_step_1);%In ms

%% PE TRAJECTORY
offset = -0;%for siemens raw data, it was -1. If 0, it does nothing
Lines = mod( single(meas.head.idx.kspace_encode_step_1)-1 + offset, NY(2) ) + 1;
Partitions = mod( single(meas.head.idx.kspace_encode_step_2) -1 +offset, NY(3) ) +1;

kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE 
rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices

%% SAVE
[pathOuTemp,rec.Names.Name] = fileparts(fileName);
if isempty(pathOu); rec.Names.pathOu = pathOuTemp;end

%%% WRITE RAW
if writeRAW
    save(strcat(fileName,'.mat'),'rec','-v7.3');    
    fprintf('File saved:\n   %s\n', fileName);
end

%%% WRITE METADATA TO JSON
recJSON = rec;
%Remove all memory intensive arrays
recJSON.y=[];
recJSON.S=[];recJSON.N=[];recJSON.Assign=[];
recJSON.x=[];recJSON.M=[];
%Save JSON
savejson('',recJSON,sprintf('%s.json',fileName));%writeJSON has specific fields 

end