

function rec = resampleRec (rec, resNew, supportReadout, underSampling)

%RESAMPLEREC  Resamples a reconstruction structure to a new resolution and under-sampling pattern. For original structure creation, see TWIX2Rec.m
%   REC=RESAMPLEREC(REC,{RESNEW},{SUPPORTREADOUT},{UNDERSAMPLING})
%   * REC is the original reconstruction structure.
%   * {RESNEW} is the resolution at which to resample the reconstruction structure.
%   * {SUPPORTREADOUT} is the support in the readout in percentage to extract. Especially useful to reconstruct a thin slice of a high-res volmetric acquisition.
%   * {UNDERSAMPLING} is the equidistant under-sampling factor in both PE dimensions (defaults to [1 1]).
%   ** REC is the modified reconstruction structure.
%

if nargin < 2  || isempty(resNew);resNew=[];end
if nargin < 3  || isempty(supportReadout);supportReadout=[];end
if nargin < 4  || isempty(underSampling);underSampling=[];end

%% PREPARE PARAMETERS
if isempty(resNew) && isempty(supportReadout) && isempty(underSampling);return;end

%%% Resolution
resOld=rec.Enc.AcqVoxelSize;
if length(resNew)==1; resNew = resNew * ones([1 3]); else; resNew(end+1:3) = resOld((length(resNew)+1):3);end

%%% Undersampling
if length(underSampling)==1; underSampling = underSampling * ones([1 2]);end

%%% Image and k-space array size
NImageOrig = rec.Enc.FOVSize(1:3);
NKSpaceOrig = size(rec.y); NKSpaceOrig=NKSpaceOrig(1:3);

wasOverSampled = rec.Enc.FOVSize==.5*rec.Enc.AcqSize;
isPT = isfield(rec,'PT');%Whether reconstruction structure contains Pilot Tone signal

%% CHANGING RESOLUTION
if ~isempty(resNew) && ~isequal(resNew, resOld)
    
    NKSpaceNew = round(NKSpaceOrig .* (resOld./resNew));
    resNew = resOld .* (NKSpaceOrig./NKSpaceNew);
    NImageNew = round(rec.Enc.FOVmm./resNew); %TODO: sort out what values to give this
    warning('TODO: sort out what values to give NImageNew.\n')
    
    %%% K-SPACE
    fprintf('Resampling k-space to a different resolution.\n');
    rec.y = resampling(rec.y, NKSpaceNew); %rec.y is in the image domain
    
    rec.Enc.AcqSize=NKSpaceNew(1:3);
    if wasOverSampled; rec.Enc.AcqSize(1) = 2*rec.Enc.AcqSize(1);end

    rec.Enc.FOVSize= NImageNew(1:3);

    %%% PILOT TONE SIGNAL 
    if isPT
        fprintf('Resampling PT signal to a different resolution.\n');
        NSlice = size(rec.PT.pSliceImage);
        NSliceNew = NSlice ; NSliceNew(2:3) = NKSpaceNew(2:3);
        for n=2:3;rec.PT.pSliceImage = fftshiftGPU(rec.PT.pSliceImage,n);end%The extraction in the fft should be in the shifted case since that is where sampling is defined in
        rec.PT.pSliceImage = resampling(rec.PT.pSliceImage, NSliceNew(1:3) );
        for n=2:3;rec.PT.pSliceImage = ifftshiftGPU(rec.PT.pSliceImage,n);end
        rec.PT.idxRO=[];
    end

    %%% PARAMETERS
    %Orientation
    fprintf('Modify NIFTI header for altered resolution.\n');
    [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize,rec.Par.Mine.APhiRec,'resampling',[],NImageOrig, NImageNew);
    rec.Par.Mine.Asca=[];%Set empty since not compatible with new resolution (but normally not needed anyway)
    rec.Par.Mine.Atra=[];
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'ipermute',rec.Par.Mine.permuteHist{1});
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'resampling',[],NImageOrig, NImageNew);
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'permute',rec.Par.Mine.permuteHist{1});

    rec.Enc.FOVmm = rec.Enc.AcqVoxelSize .* rec.Enc.FOVSize; 
    
    %Sampling trajectory
    fprintf('Modify sampling parameters for altered resolution.\n');
    z = rec.Assign.z;

        %Shift to make compatible with DISORDER processing
        kShift = floor((NKSpaceOrig)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
        z{2}= z{2} + kShift(2); % 2nd PE 
        z{3}= z{3} + kShift(3); % 3rd PE = slices

        %Find centre and what samples to extract
        zeroF=ceil((NKSpaceOrig+1)/2);%Nmin=min(NOld,NNew);Nmax=max(NOld,NNew);
        orig=zeroF-ceil((NKSpaceNew-1)/2); %YB: These are the elements to extract from the bigger F matrix to resample the array
        fina=zeroF+floor((NKSpaceNew-1)/2);

        idxSamples = (z{2} >= orig(2)) & (z{2} <= fina(2)) & (z{3} >= orig(3)) & (z{3} <= fina(3));

        %Extract samples
        z{2} = dynInd( z{2} , idxSamples, 2);
        z{3} = dynInd( z{3} , idxSamples, 2);

        %Add offset to samples to make first index=1
        z{2} = z{2} - (orig(2)-1);
        z{3} = z{3} - (orig(3)-1);

        %Add shift for DISORDER again:
        kShift = floor((NKSpaceNew)/2) + 1;
        z{2} = z{2} - kShift(2);
        z{3} = z{3} - kShift(3);

        %Re-assign
        rec.Assign.z = z; z = [];
    
    %Sequence parameters
    fprintf('Modify sequence parameters for altered resolution.\n');
    rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / NImageNew(1);%In Hz/voxel in the readout direction
    rec.Par.Labels.ScanDuration=inf;%This would not be reliable information.
    rec.Par.Labels.TFEfactor = length(rec.Assign.z{2});
    NY = NKSpaceNew(1:3);
    rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};

    if isfield(rec.Enc,'DISORDER'); rec.Enc.DISORDER.Warning='Parameters deprecated after k-space resampling.';end
end

%% CHANGING READOUT SUPPORT
NY = size(rec.y);
if ~isempty(supportReadout)
    fprintf('Extracting part of FOV with range [%.2f --> %.2f].\n', supportReadout(1), supportReadout(2));
    %Extraction
    vr = max(1,round(NY(1)*supportReadout(1))) : min(NY(1),round(NY(1)*supportReadout(2)));
    idx = {vr,':',':',':'};
    rec.y = dynInd(rec.y,idx ,1:4);
    %Geomtry
    [~,rec.Par.Mine.APhiRec]=mapNIIGeom([],rec.Par.Mine.APhiRec,'dynInd',idx,NY,size(rec.y));    
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'ipermute',rec.Par.Mine.permuteHist{1});
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'dynInd',idx,NY,size(rec.y));    
    [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'permute',rec.Par.Mine.permuteHist{1});
    %Labels
    rec.Par.Labels.Bandwidth = rec.Par.Labels.voxBandwidth * size(rec.y,1);%In Hz/voxel in the readout direction
    %Array sizes
    NY=size(rec.y);
    rec.Enc.AcqSize=NY(1:3);
    if wasOverSampled; rec.Enc.AcqSize(1) = 2*rec.Enc.AcqSize(1);end
    
    rec.Enc.FOVSize(1)=NY(1);
    if isfield(rec.Enc,'FOVmm'); rec.Enc.FOVmm = rec.Enc.FOVSize.* rec.Enc.AcqVoxSiz; end
end

%% CHANGE UNDER-SAMPLING
if ~isempty(underSampling) && ~isequal(underSampling,[1 1])
    
    R = underSampling;%Easier to read
    assert(all(mod(R,1)==0), 'resampleRec:: Undersampling only allowed to be integers.')
    
    %%% Changing rec.Assign.z
    fprintf('Modify sampling parameters for altered undersampling.\n');
    z = rec.Assign.z;
    [kIndex, idxToSort, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    NY = size(rec.y);
    %%% Shift to make compatible with DISORDER recon
    kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
    z{2}= z{2} + kShift(2); % 2nd PE 
    z{3}= z{3} + kShift(3); % 3rd PE = slices
    
    %Select indices that are 
    idxToExtract = cell(1,2); %PE dimensions
    startIdx = [1 1];%Max value = R
    for PEDim=1:2; idxToExtract{PEDim} = startIdx(PEDim):R(PEDim):NY(PEDim+1);end%+1 since first dim is the RO
    
    gridToExtract = zeros(NY(2:3),'like',real(rec.y));
    gridToExtract = dynInd(gridToExtract, idxToExtract ,1:2, 1);%Set value to 1 that are extracted
    
    for n=2:3
        %rec.y=fftshiftGPU(rec.y,n);
        rec.y=fftGPU(rec.y,n)/NY(n);%Divide by NY to make array size not affect PT signal
        rec.y=fftshiftGPU(rec.y,n);%fftshift since in solveXT there is an ifft on timeIndex
        
        if isPT;rec.PT.pSliceImage = fftshiftGPU(rec.PT.pSliceImage,n);end
        if isPT;rec.PT.pSliceImage = fftGPU(rec.PT.pSliceImage,n)/NY(n);end
        if isPT;rec.PT.pSliceImage = fftshiftGPU(rec.PT.pSliceImage,n);end
    end
    fprintf('Modify k-space for altered undersampling.\n');
    removeSamples = 1;
    if removeSamples; rec.y = dynInd(rec.y,idxToExtract,2:3);end
    NYNew = size(rec.y);
    if isPT
        fprintf('Modify Pilot Tone signal for altered undersampling.\n');
        if removeSamples; rec.PT.pSliceImage = dynInd(rec.PT.pSliceImage,idxToExtract,2:3);end
    end
    
    rec.Enc.AcqSize=NYNew(1:3);
    if wasOverSampled; rec.Enc.AcqSize(1) = 2*rec.Enc.AcqSize(1);end
    
     for n=2:3
        rec.y=ifftshiftGPU(rec.y,n);
        rec.y=ifftGPU(rec.y,n)*(NYNew(n));
        %rec.y=ifftshiftGPU(rec.y,n);
        
        if isPT;rec.PT.pSliceImage = ifftshiftGPU(rec.PT.pSliceImage,n);end
        if isPT;rec.PT.pSliceImage = ifftGPU(rec.PT.pSliceImage,n)*NYNew(n);end
        if isPT;rec.PT.pSliceImage = ifftshiftGPU(rec.PT.pSliceImage,n);end
     end
    
    
    gridToExtract = resSub(gridToExtract, 1:2);
    idxSamples = gridToExtract==1;

    idxSamples = squeeze(dynInd(idxSamples, idxToSort, 1)).';%sorted
    
    %%% Extract samples
    z{2} = dynInd( z{2} , idxSamples, 2);
    z{3} = dynInd( z{3} , idxSamples, 2);

    if removeSamples
        %%% Make jumps in k-space coordinates 1 instead of R
        z{2} = ceil( z{2}/R(1) ); %No +1 since if startIdx = R, you don't want first index to be 2
        z{3} = ceil( z{3}/R(2) ); 
    end
     
    %%% Add shift for DISORDER again with update array size
    kShift = floor((NYNew)/2) + 1;
    z{2} = z{2} - kShift(2);
    z{3} = z{3} - kShift(3);

    %%% Re-assign
    rec.Assign.z = z;z = [];

    rec.Par.Labels.TFEfactor = length(rec.Assign.z{2});
    rec.Enc.kRange={[-NYNew(1) NYNew(1)-1],[-NYNew(2)/2 NYNew(2)/2-1],[-NYNew(3)/2 NYNew(3)/2-1]}; 
    
    if isfield(rec.Enc,'DISORDER'); rec.Enc.DISORDERInfo.Warning='Parameters deprecated after k-space resampling.';end
end

%% (OPTOINAL) RE-SORT PT SIGNAL
if isPT
    tt = rec.PT.pSliceImage;
    NY = size(rec.y);
    for n=2:3
        tt=fftshiftGPU(tt,n);
        tt=fftGPU(tt,n)/NY(n);%Divide by NY to make array size not affect PT signal
        tt=fftshiftGPU(tt,n);%fftshift since in solveXT there is an ifft on timeIndex
    end

    [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);

    tt = resSub(tt, 2:3);
    rec.PT.pTimeTest = squeeze(dynInd(tt, idx, 2)).';%sorted

    coilPlot=1:min(10,size(rec.PT.pTimeTest,1));
    visPTSignal (dynInd(rec.PT.pTimeTest,coilPlot,1), [], [],[], 2, [], [], 201, replace(sprintf('\\textbf{Temporal PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [\#TR]');    

    h=figure('color','w'); imshow(timeMat , []);title('Time indices of k-space sampling')
    set(h,'color','w','Position',get(0,'ScreenSize'));

end