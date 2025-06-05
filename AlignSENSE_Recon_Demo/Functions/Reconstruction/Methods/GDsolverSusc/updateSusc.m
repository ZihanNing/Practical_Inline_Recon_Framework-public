
function E = updateSusc (E, x, N, MS, MT, parXB, Labels)

if nargin<4 || isempty(MS);MS=ones(1,3);end
if nargin<5 || isempty(MT);MT=eye(4);end
isValidGeom(MS,MT,2);
gpuIn=single(gpuDeviceCount);

if length(MS)==1; MS = MS*ones([1,length(N)]);end

if ~isfield(E, 'Ds')
    fprintf('updateSusc:: Creating susceptibility coeffients for B0 modelling\n')
    
    % Initialise
    yawIdx = 4; rollIdx=5; pitchIdx=6;%Reverse order as in sincRigidTransformation convention/parameters in RAS
    [E.Ds.f, E.Ds.fDeriv] = higherOrderFunHandle(1 , yawIdx, rollIdx , pitchIdx);
    E.Ds.chi =  zeros([N , 1, 1, 1], 'like',real(x) );
    E.Ds.padDim_mm = 0*[20 20 20];%In mm
        E.Ds.padDim = round(E.Ds.padDim_mm./MS);%From mm to voxels
    E.Ds.chi = padArrayND( E.Ds.chi , E.Ds.padDim, 1, 0, 'both');
    Nchi = size(E.Ds.chi);
    Nkernels = N + 2*E.Ds.padDim;
    E.Ds.Nchi = Nchi;
    
    % Parameters
    E.Ds.TE = Labels.TE * 1e-3;%in seconds  
    if ~isfield( Labels , 'H0'); E.Ds.H0 = 7; fprintf('updateSusc:: H0 set to 7T.\n');else E.Ds.H0=Labels.H0;end %Tesla
    if ~isfield( Labels , 'gamma'); E.Ds.gamma = 42.58; fprintf('updateSusc:: gamma set to 42.58 Mhz/T.\n');else E.Ds.gamma=Labels.gamma;end %MHz/T

    if isfield(parXB,'chiKernel') && isfield(parXB.chiKernel,'filter')  && isfield(parXB.chiKernel.filter,'sp');E.Ds.kernelFilter.sp = .8; fprintf('updateSusc:: kernelFiltSp set to 0.8.\n'); else E.Ds.kernelFilter.sp = parXB.chiKernel.filter.sp; end
    if isfield(parXB,'chiKernel') && isfield(parXB.chiKernel,'filter')  && isfield(parXB.chiKernel.filter,'gibbs');E.Ds.kernelFilter.gibbs = .3; fprintf('updateSusc:: kernelFiltGibbs set to 0.3.\n'); else E.Ds.kernelFilter.gibbs = parXB.chiKernel.filter.gibbs; end

    if isfield(parXB,'orderSusc'); E.Ds.orderSusc = parXB.orderSusc; else E.Ds.orderSusc=[0 1];end
    
    % Regularisation
    E.Ds.Reg = parXB.Reg;
    
    % Optimisation
    E.Ds.Optim = parXB.Optim;
    if ~isfield( E.Ds.Optim, 'alphaList'); E.Ds.Optim.alphaList=10.^(-1* [-2:0.7:6]);end
    if strcmp(Labels.FatShiftDir,'F');E.Ds.Optim.Mbg=1:floor((1-parXB.redFOV)*N(3));elseif strcmp(Labels.FatShiftDir,'H');E.Ds.Optim.Mbg=floor(1+parXB.redFOV*N(3)):N(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Ds.Optim.Mbg=1:N(3);end
    
else %Need to resample the susceptibility distribution
    fprintf('updateSusc:: Updating susceptibility coeffients for B0 modelling\n')
    
    E.Ds.padDim = round(E.Ds.padDim_mm./MS);%From mm to voxels
    Nchi = N + 2*E.Ds.padDim;
    Nkernels = N + 2*E.Ds.padDim;
    %%% Resample
    E.Ds.chi=resampling(E.Ds.chi,Nchi,0,2*ones(1,3));  %Resample in DCT domain      
end

% Constrain  -MOVE TO ABOVE
if isfield(parXB ,'filterChi')
    if isfield ( parXB.filterChi, 'sp'); E.Ds.C.filterParams.sp = parXB.filterChi.sp; else E.Ds.C.filterParams.sp = 10;end
    if isfield ( parXB.filterChi, 'gibbsRinging'); E.Ds.C.filterParams.gibbsRinging = parXB.filterChi.gibbsRinging; else E.Ds.C.filterParams.gibbsRinging = 0;end
    E.Ds.C.H = buildFilter(2*Nchi(1:3),'tukeyIso',(MS/E.Ds.C.filterParams.sp)*[1 1 1],gpuIn,E.Ds.C.filterParams.gibbsRinging,1);%Filter in Cosine domain 
end
    
if strcmp(Labels.FatShiftDir,'F');E.Ds.Optim.Mbg=1:floor((1-parXB.redFOV)*N(3));elseif strcmp(Labels.FatShiftDir,'H');E.Ds.Optim.Mbg=floor(1+parXB.redFOV*N(3)):N(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.Ds.Optim.Mbg=1:N(3);end

%%% CREATE K-SPACE COORDINATES
FOV = Nkernels(1:3).*MS;
type = class(gather(x));
Nkernels=cast(Nkernels, type);
FOV=cast(FOV, type);
[kx, ky, kz]=linearTerms(Nkernels(1:3), Nkernels(1:3)/2./FOV);%*N/2 to compensate for scaling in linearTerms.m and ./FOV to account for anisotropy

%%% TAKE ORIENTATION INTO ACCOUNT 
[kx, ky, kz] = rotateCoordinates(kx, ky, kz, inv(MT));

%%% MAKE KERNELS
%Zero order term
if ismember(0,E.Ds.orderSusc)
    E.Ds.kernelStruct.Kequilibrium = 1/3 - kz.^2./(kx.^2+ky.^2+kz.^2);
    E.Ds.kernelStruct.Kequilibrium(isnan(E.Ds.kernelStruct.Kequilibrium))=1/3;%If 0, mean induced field=0, if 1/3, mean induced field = 1/3*mean(susc)
else
    E.Ds.kernelStruct.Kequilibrium=[];
end

%First order term
if ismember(1,E.Ds.orderSusc)
    E.Ds.kernelStruct.Kroll = -2*kx.*kz./(kx.^2+ky.^2+kz.^2);   E.Ds.kernelStruct.Kroll(isnan(E.Ds.kernelStruct.Kroll))=0;
    E.Ds.kernelStruct.Kpitch = 2*ky.*kz./(kx.^2+ky.^2+kz.^2);   E.Ds.kernelStruct.Kpitch(isnan(E.Ds.kernelStruct.Kpitch))=0;
else
    E.Ds.kernelStruct.Kroll=[];
    E.Ds.kernelStruct.Kpitch=[];
end

%E.Ds.kernelStruct.Kequilibrium = ones(size(kz),'like', kz);
 
%%% FILTERS
%Kernel 
E.Ds.kernelStruct.H = buildFilter( Nkernels(1:3),'tukeyIso',E.Ds.kernelFilter.sp,gpuIn,E.Ds.kernelFilter.gibbs );
E.Ds.kernelStruct.H = ones(size(kz),'like', kz);

% Chi constraints
% TO DO

%%% CREATE LCMs and B0
E = convertChiToTaylorMaps(E);


