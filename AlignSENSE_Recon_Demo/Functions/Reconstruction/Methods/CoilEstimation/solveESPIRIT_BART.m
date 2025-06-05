
function rec = solveESPIRIT_BART(rec)

if ~isfield(rec.Alg.parE,'BARTPath'); rec.Alg.parE.BARTPath='/home/ybr19/Software/gitRepos/BART/bart-0.7.00';end

%%% Add BART to the path
addpath(genpath(rec.Alg.parE.BARTPath))

%%% Define constansts and variables
parE = rec.Alg.parE;
N = size(rec.y);
ND=3;%Number of dimensions for estimation

%%% Define calibration region
DeltaK=(1./(N(1:ND).*rec.Enc.AcqVoxelSize(1:ND)));%1/FOV - YB: in units of 1/mm
DeltaK=DeltaK(1:ND);

KMaxCalib = 1./(2*rec.Alg.parE.NC); %YB: in units of 1/mm
parE.NC=2*KMaxCalib./DeltaK;%Resolution (array size) of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution %Now parE.NC is in units of voxels
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
%YB:Now NC is in units of voxels to include in AC region and not resolution anymore

%%% Define kernel size
parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles 
parE.K=max(parE.K,parE.Kmin);
parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution
    
 
%%% Call ESPIRiT from BART toolbox
num_acs = parE.NC;       % size of calibration k-space fopr ESPIRiT
k = parE.K;              % kernel size
%c = 1;                   % threshold for sensitivity mask    

dataESPIRiT = bart('fft 7', single( cat(4, rec.x, rec.y) ));%go to k-space for ecalib funcitonality
[coilSensitivities, eigenMaps] = bart( ['ecalib -r ', num2str(num_acs),' -k ', num2str(k)], dataESPIRiT);

%%% Store in the structure
rec.S = dynInd(coilSensitivities, 1, 5);%First eigenmap
rec.S = dynInd(rec.S, 2:size(rec.S,4), 4);%Exclude body coil

rec.W = eigenMaps;
rmpath(genpath(rec.Alg.parE.BARTPath))

