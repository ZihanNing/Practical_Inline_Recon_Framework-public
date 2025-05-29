function rec=reconInvert_gadg(rec,typ)

% ZN: this is a gadgetron version of the Invert
% modified based on Lucilio's reconInvert (which can be found in
% Functions_SENSE) to adapt the input of ISMRMRD dataset
% Work with Gadgetron to achieve inline reconstruction
%
% 27-Jan-25
% Zihan Ning

%RECONINVERT   Inverts data
%   REC=RECONINVERT(REC,TYP,{FIELD})
%   * REC is a recon structure
%   * TYP is the type of data to invert, one of the following, 1->body, 
%2->surface, 3->data
%   ** REC is a recon structure

typV=['B' 'S' 'x'];
t=typV(typ);
% TW=rec.TW;%This is the data to invert
if typ<=2;supportReadout=rec.Alg.supportReadoutS;else supportReadout=rec.Alg.supportReadout;end

ND=16;
gpu=useGPU;

hdrUPDouble_l = length(rec.ismrmrd.hdr.hdrUPDouble);
hdrUPLong_l = length(rec.ismrmrd.hdr.hdrUPLong);
hdrUPString_l = length(rec.ismrmrd.hdr.hdrUPString);

%ACQUIRED GRID SIZES
switch t
    case 'x'
        if length(size(rec.ismrmrd.data))>=4; rec.Enc.(t).ThreeD = 1; end % ZN: whether 3D imaging
    case 'S'
        if length(size(rec.ismrmrd.refdata))>=4; rec.Enc.(t).ThreeD = 1; end % ZN: whether 3D imaging
end
% NImageCols=TW.hdr.Config.NImageCols;%Alternative would be TW.hdr.Meas.NImageCols
% if isfield(TW.hdr.Meas,'NImageCols');NImageCols=TW.hdr.Meas.NImageCols;else NImageCols=TW.(field).NCol;end
NImageCols = rec.ismrmrd.hdr.hdrConnection.encoding.reconSpace.matrixSize.x;
% if TW.(field).NPar==1;rec.Enc.(t).AcqN=[NImageCols TW.(field).NLin-TW.(field).skipLin TW.(field).NSli];rec.Enc.(t).ThreeD=0;%MS2D
% else rec.Enc.(t).AcqN=[NImageCols TW.(field).NLin-TW.(field).skipLin TW.(field).NPar-TW.(field).skipPar];rec.Enc.(t).ThreeD=1;%3D
% end
% if rec.Enc.(t).ThreeD;fprintf('Volumetric encoding\n');else fprintf('Multi-slice encoding\n');end
switch t
    case 'x'
        rec.Enc.(t).AcqN = [2 * NImageCols,size(rec.ismrmrd.data,2) size(rec.ismrmrd.data,3)]; % datasize after zero-padding (due to PF for example)
    case 'S'
        rec.Enc.(t).AcqN = [2 * NImageCols size(rec.ismrmrd.refdata,2) size(rec.ismrmrd.refdata,3)];
end
fprintf('Acquired size (no PF):%s\n',sprintf(' %d',rec.Enc.(t).AcqN));

%NUMBER OF COILS
% rec.Enc.(t).CoilN=TW.(field).NCha;
switch t
    case 'S' % refscan
        rec.Enc.(t).CoilN = size(rec.ismrmrd.refdata,4);
    case 'x' % high-res
        rec.Enc.(t).CoilN = size(rec.ismrmrd.data,4);
end
fprintf('Number of coils: %s\n',sprintf('%d',rec.Enc.(t).CoilN));

%OVERSAMPLING FACTORS
rec.Enc.(t).Overs=ones(1,3);%Readout is treated indepedently in function invert readout
% if ~isempty(TW.hdr.Meas.flPhaseOS);rec.Enc.(t).Overs(2)=1+TW.hdr.Meas.flPhaseOS;end
% if ~isempty(TW.hdr.Meas.flSliceOS);rec.Enc.(t).Overs(3)=1+TW.hdr.Meas.flSliceOS;end
for i=1:hdrUPDouble_l
    if strcmp(rec.ismrmrd.hdr.hdrUPDouble(i).name,'flPhaseOS_zihan')
        rec.Enc.(t).Overs(2) = 1+rec.ismrmrd.hdr.hdrUPDouble(i).value;
    end
end
for i=1:hdrUPDouble_l
    if strcmp(rec.ismrmrd.hdr.hdrUPDouble(i).name,'flSliceOS_zihan')
        rec.Enc.(t).Overs(3) = 1+rec.ismrmrd.hdr.hdrUPDouble(i).value;
    end
end
fprintf('Oversampling factors:%s\n',sprintf(' %.3f',rec.Enc.(t).Overs));

%ACCELERATION FACTORS
rec.Enc.(t).Accel=ones(1,3);%Readout is treated indepedently
% rec.Enc.(t).RegridMode=TW.hdr.Meas.alRegridMode(1);
rec.Enc.(t).RegridMode = 1; % ZN: double check whether it is suitable for highres as well
fprintf('Regrid mode: %s\n',sprintf('%d',rec.Enc.(t).RegridMode));
% rec.Enc.(t).Accel(2)=1/TW.hdr.MeasYaps.sKSpace.dPhaseResolution;
% rec.Enc.(t).Accel(3)=1/TW.hdr.MeasYaps.sKSpace.dSliceResolution;
for i=1:hdrUPDouble_l
    if strcmp(rec.ismrmrd.hdr.hdrUPDouble(i).name,'dPhaseResolution_zihan')
        rec.Enc.(t).Accel(2) = 1/rec.ismrmrd.hdr.hdrUPDouble(i).value;
    end
end
for i=1:hdrUPDouble_l
    if strcmp(rec.ismrmrd.hdr.hdrUPDouble(i).name,'dSliceResolution_zihan')
        rec.Enc.(t).Accel(3) = 1/rec.ismrmrd.hdr.hdrUPDouble(i).value;
    end
end
fprintf('Acceleration factors:%s\n',sprintf(' %.3f',rec.Enc.(t).Accel));
% if rec.Enc.(t).RegridMode<=1;rec.Enc.(t).AcqN(1)=rec.Enc.(t).AcqN(1)/TW.hdr.Dicom.flReadoutOSFactor;end
for i=1:hdrUPDouble_l
    if strcmp(rec.ismrmrd.hdr.hdrUPDouble(i).name,'flReadoutOSFactor_zihan')
        flReadoutOSFactor = rec.ismrmrd.hdr.hdrUPDouble(i).value;
    end
end
if rec.Enc.(t).RegridMode<=1;rec.Enc.(t).AcqN(1)=rec.Enc.(t).AcqN(1)/flReadoutOSFactor;end

%DEAL WITH BIPOLAR READOUT (ONLY FOR GADGETRON)
% ZN: I suppose that this will be automatic handled if you use mapVBVD &
% twix
if typ == 3
    flag_bipolar = 1; % ZN: good to replace this with a ismrmrd header
    if flag_bipolar
        data = rec.ismrmrd.data;
        Necho = size(data,ndims(data));
        if Necho>=2;evenEchoIdx = 2:2:Necho;end % ZN: only flip the even echo
        data(:,:,:,:,evenEchoIdx) = flip(data(:,:,:,:,evenEchoIdx),1);
        rec.ismrmrd.data = data; 
    end
end

%PARTIAL FOURIER
% ZN: the pre-processing gadget already did the zero-padding (due to PF),
% therefore no need to do it here
if typ == 3 % ZN: if high-res scan (only need to calc APF, the logic matrix of padded area, for mask generation) 
%     rec.Enc.(t).AcqNNoPF=rec.Enc.(t).AcqN;
%     if rec.Enc.(t).ThreeD;rec.Enc.(t).AcqN=[rec.Enc.(t).AcqN(1) TW.hdr.Config.N0FPEFTLength TW.hdr.Config.N0F3DFTLength];
%     else rec.Enc.(t).AcqN=[rec.Enc.(t).AcqN(1) TW.hdr.Config.N0FPEFTLength rec.Enc.(t).AcqN(3)];
%     end
%     pad=rec.Enc.(t).AcqN-rec.Enc.(t).AcqNNoPF;
%     rec.Enc.(t).APF=cell(1,3);
%     for d=1:3
%         rec.Enc.(t).APF{d}=ones(rec.Enc.(t).AcqNNoPF(d),1,'single');
%         if useGPU;rec.Enc.(t).APF{d}=gpuArray(rec.Enc.(t).APF{d});end
%         rec.Enc.(t).APF{d}=padarray(rec.Enc.(t).APF{d},pad(d),0,'post');
%     end   
%     pad=[pad(1) 0 pad(2:3)];
%     rec.TWD.(t)=padarray(rec.TWD.(t),pad,0,'post');
%     if rec.Enc.(t).ThreeD;dims=[3 4];else dims=3;end
%     cshift=zeros(1,4);
%     for d=dims
%         if d==3;cshift(3)=ceil((rec.Enc.(t).AcqN(2)+1)/2)-TW.(field).centerLin(1);
%         else cshift(4)=ceil((rec.Enc.(t).AcqN(3)+1)/2)-TW.(field).centerPar(1);
%         end                
%         rec.Enc.(t).APF{d-1}=circshift(rec.Enc.(t).APF{d-1},cshift(d));
%     end
%     rec.TWD.(t)=circshift(rec.TWD.(t),cshift);
%     for d=1:3;rec.Enc.(t).APF{d}=shiftdim(rec.Enc.(t).APF{d},d-1);end
    rec.Enc.(t).APF = findMatrixByNonZero(rec.ismrmrd.data);
end

%RECONSTRUCTED GRID SIZES
rec.Enc.(t).RecN=round(rec.Enc.(t).AcqN.*rec.Enc.(t).Accel);%Oversampling will be removed at the end
fprintf('Reconstructed size:%s\n',sprintf(' %d',rec.Enc.(t).RecN));

%OUTPUT GRID SIZES
rec.Enc.(t).OutN=round(rec.Enc.(t).AcqN.*rec.Enc.(t).Accel./rec.Enc.(t).Overs);%Oversamling will be removed at the end
fprintf('Output size:%s\n',sprintf(' %d',rec.Enc.(t).OutN));%This is not used

%SLICE SEPARATION
% aux=TW.hdr.MeasYaps.sSliceArray.asSlice{1};
rec.Enc.(t).SliceSeparation=0;
if 0 %~rec.Enc.(t).ThreeD % skipped for both high-res & refscan
    posSl1=[aux.sPosition.dSag aux.sPosition.dCor aux.sPosition.dTra];
    aux=TW.hdr.MeasYaps.sSliceArray.asSlice{2};
    posSl2=[aux.sPosition.dSag aux.sPosition.dCor aux.sPosition.dTra];
    rec.Enc.(t).SliceSeparation=norm(posSl1-posSl2);
    fprintf('Slice separation: %s\n',sprintf('%.2f',rec.Enc.(t).SliceSeparation));
end


%FULL ACQUIRED FOV
% aux=TW.hdr.MeasYaps.sSliceArray.asSlice{1};
% if rec.Enc.(t).ThreeD
%     rec.Enc.(t).AcqFOV=[aux.dReadoutFOV aux.dPhaseFOV aux.dThickness];rec.Enc.(t).SliceThickness=0;
% else 
%     rec.Enc.(t).AcqFOV=[aux.dReadoutFOV aux.dPhaseFOV rec.Enc.(t).SliceSeparation*rec.Enc.(t).AcqN(3)];rec.Enc.(t).SliceThickness=aux.dThickness;
%     fprintf('Slice thickness: %.2f\n',rec.Enc.(t).SliceSeparation);
% end
dReadoutFOV = rec.ismrmrd.hdr.hdrConnection.encoding.reconSpace.fieldOfView_mm.x;
dPhaseFOV = rec.ismrmrd.hdr.hdrConnection.encoding.reconSpace.fieldOfView_mm.y;
dThickness = rec.ismrmrd.hdr.hdrConnection.encoding.reconSpace.fieldOfView_mm.z;
rec.Enc.(t).AcqFOV = [dReadoutFOV , dPhaseFOV , dThickness];
rec.Enc.(t).AcqFOV=rec.Enc.(t).AcqFOV.*rec.Enc.(t).Overs;
fprintf('Acquired field of view:%s\n',sprintf(' %.2f',rec.Enc.(t).AcqFOV));

%ACQUIRED VOXEL SIZE
rec.Enc.(t).AcqDelta=rec.Enc.(t).AcqFOV./rec.Enc.(t).AcqN(1:3);
fprintf('Sampled resolution:%s\n',sprintf(' %.2f',rec.Enc.(t).AcqDelta));

%TO EXTRACT FOV
if ~isempty(supportReadout);vr=round(rec.Enc.(t).AcqN(1)*supportReadout(1)):round(rec.Enc.(t).AcqN(1)*supportReadout(2));else vr=1:rec.Enc.(t).AcqN(1);end
%if ~isempty(supportReadout);vr=round(rec.Enc.(t).RecN(1)*supportReadout(1)):round(rec.Enc.(t).RecN(1)*supportReadout(2));else vr=1:rec.Enc.(t).RecN(1);end

%SLICE AND PHASE ORDER
switch t
    case 'S'
        rec.TWD.(t) = rec.ismrmrd.refdata; rec.ismrmrd.refdata = [];
    case 'x'
        rec.TWD.(t) = rec.ismrmrd.data; rec.ismrmrd.data = [];       
end
% ZN: be careful, the next step might just work for the order of MEGE
% currently
rec.TWD.(t) = permute(rec.TWD.(t),[1,4,2,3,6,7,8,5]); % ZN: to make the dimension of MRD to match the dimension of twix

NX=size(rec.TWD.(t));

if 0 %~rec.Enc.(t).ThreeD
    [rec.Enc.(t).SliceOrder,rec.Enc.(t).PEOrder,rec.Enc.(t).ShotOrder,rec.Enc.(t).EchoOrder,rec.Enc.(t).SlicePositions,rec.Enc.(t).SliceOffsets,iSlicePositionsMin,rec.Enc.(t).SliceOrderReduced,rec.Enc.(t).ShotOrderReduced]=computeSlicePEShotAndEchoOrders(TW.(field),NX);
else
    iSlicePositionsMin=1;
end
if ~isempty(rec.Alg.maxNumberRepeats);rec.TWD.(t)=dynInd(rec.TWD.(t),1:min(size(rec.TWD.(t),9),rec.Alg.maxNumberRepeats),9);end

%REGRIDDING 
% ZN: this is skipped for both high-res and refscan
if rec.Enc.(t).RegridMode>1%From https://github.com/wtclarke/pymapvbvd/blob/master/mapvbvd/read_twix_hdr.py line 206
    ncol=TW.hdr.Meas.alRegridDestSamples(1);%200
    dwelltime=TW.hdr.Meas.aflRegridADCDuration(1)/ncol;%560/200
    start=TW.hdr.Meas.alRegridDelaySamplesTime(1);%60
    rampup_time=TW.hdr.Meas.alRegridRampupTime(1);%190
    flattop_time=TW.hdr.Meas.alRegridFlattopTime(1);%300
    rampdown_time=TW.hdr.Meas.alRegridRampdownTime(1);%190

%     fprintf('Number of cols: %d\n',ncol);
%     fprintf('Dwell time: %.3f\n',dwelltime);
%     fprintf('Start: %.3f\n',start);
%     fprintf('Rampup time: %.3f\n',rampup_time);
%     fprintf('Flattop time: %.3f\n',flattop_time);
%     fprintf('Rampdown time: %.3f\n',rampdown_time);
%BandwidthPerPixelPhaseEncode
%EchoSpacing_eff_us
%EchoSpacing_us
%MaxwellIntegralROGradient
%alTE
%flBandwidthPerPixelPhaseEncode
%lEchoSpacing

    gr_adc=zeros([ncol 1],'single');
    time_adc=start+dwelltime*(0.5+0:ncol-1)';
    ixUp=time_adc<rampup_time;
    ixFlat=time_adc<=rampup_time+flattop_time & ~ixUp;
    ixDn=~ixUp & ~ixFlat;
    gr_adc(ixFlat)=1;
    if rec.Enc.(t).RegridMode==2%Trapezoidal gradient
        gr_adc(ixUp)=time_adc(ixUp)/rampup_time;
        gr_adc(ixDn)=1-(time_adc(ixDn)-rampup_time-flattop_time)/rampdown_time;
    elseif rec.Enc.(t).RegridMode==4%Sinusoidal gradient
        gr_adc(ixUp)=sin((pi/2)*(time_adc(ixUp)/rampup_time));
        gr_adc(ixDn)=sin((pi/2)*(1+(time_adc(ixDn)-rampup_time-flattop_time)/rampdown_time));
    else
        error('Unknown regridding mode');
    end
    gr_adc=max(gr_adc,1e-4);%Make sure always positive
    rstraj=(cumtrapz(gr_adc)-ncol/2)/sum(gr_adc);
    %rec.Enc.(t).Grid=rstraj-mean(rstraj(floor(ncol/2)-1:floor(ncol/2)+1));    
    %rec.Enc.(t).Grid=(rstraj-mean(rstraj)+0.5)*ncol+1;    
    rec.Enc.(t).Grid=(rstraj-min(rstraj))/(max(rstraj)-min(rstraj))*(ncol-1)+1;    
    [A,Aw]=sincKernel(rec.Enc.(t).Grid(:),(1:ncol)',useGPU);
    rec.Enc.(t).GridMatrix=A./Aw;
else
    rec.Enc.(t).Grid=[];rec.Enc.(t).GridMatrix=[];
end

%%%%%%%%%%%%%%%%%%%%%%%GEOMETRY (based on twix)%%%%%%%%%%%%%%%%%%%%%%%
% N=rec.Enc.(t).AcqN;
% MS=rec.Enc.(t).AcqDelta;
% S=diag([MS 1]);
% rec.Geom.(t).MS=MS;
% slicePos=TW.(field).slicePos(:,iSlicePositionsMin);
% quaternionRaw = slicePos(4:7);
% R=eye(4);
% R(1:3,1:3)=convertNIIGeom(ones([1,3]),quaternionRaw','qForm','sForm');%PE-RO-SL to PCS
% R=R(:,[2 1 3 4]);
% T=eye(4);T(1:3,4)=slicePos(1:3);%in mm for center FOV (I think)
% %MT=T*R(:,[2 1 3 4])*S;%Of the center of the FOV
% MT=T*R*S;%Of the center of the FOV
% Nsub=N/2;
% if ~rec.Enc.(t).ThreeD;Nsub(3)=1/2;end
% %Nsub(1)=Nsub(1)-2;%---it is the solution
% T(:,4)=MT*[-floor(Nsub)';1];
% %T(:,4)=MT*[-Nsub';1];
% %MT=T*R(:,[2 1 3 4])*S;%Of the origin
% MT=T*R*S;
% rec.Geom.(t).patientPosition=TW.hdr.Dicom.tPatientPosition;
% fprintf('Patient position: %s\n',rec.Geom.(t).patientPosition);
% rec.Geom.(t).PCS2RAS=getPCS2RAS(rec.Geom.(t).patientPosition);
% rec.Geom.(t).MT=rec.Geom.(t).PCS2RAS*MT;
% tI=rec.Geom.(t).MT*[vr(1)-1;(rec.Enc.(t).OutN(2:3)-rec.Enc.(t).RecN(2:3))';1]; % ZN: debug, should be here, but comment out
% rec.Geom.(t).MT(1:3,4)=tI(1:3);

%%%%%%%%%%%%%%%%%%%%%%%GEOMETRY (based on ismrmrd)%%%%%%%%%%%%%%%%%%%%%%%
% ZN: calc the geom of refscan & high-res all together at the reference
% part
if typ == 2 % refscan
    % ZN: get the AcqN & AcqDelta information of high-res
    rec.Enc.x.AcqN = [2 * NImageCols,size(rec.ismrmrd.data,2) size(rec.ismrmrd.data,3)]; % Assume the NImageCols for refscan & high-res is the same
    rec.Enc.x.AcqN(1)=rec.Enc.x.AcqN(1)/flReadoutOSFactor;
    rec.Enc.x.AcqDelta=rec.Enc.S.AcqFOV./rec.Enc.x.AcqN(1:3); % the AcqFOV of refscan & high-res should be the same
    
    % ZN: to calc the geom, first permute the related parameters to [Lin RO
    % Par] first, then permute back to [RO Lin Par] after the geom calc
    rec.Enc.x.AcqDelta([1 2 3]) = rec.Enc.x.AcqDelta([2 1 3]);
    rec.Enc.x.AcqN([1 2 3]) = rec.Enc.x.AcqN([2 1 3]);

    % ZN: the permution might be tricky double check
    timeStamp = (squeeze(rec.ismrmrd.hdr.buffer.acquisition_time_stamp) );
    [timeStampSorted,idx] = sort(timeStamp(:), 'ascend');
    idx (timeStampSorted==0)=[];
    [sub1,sub2] = ind2sub(size(timeStamp), idx');
    kspace_encode_step_1 = sub1;
    kspace_encode_step_2 = sub2;
    %     Lin_ori = TW.image.Lin;
    Lin = double(kspace_encode_step_1);
    %     Par_ori = TW.image.Par;
    Par = double(kspace_encode_step_2);

    % for MEGE
    if length(NX)>=8 % ZN: might only work for MEGE dataset. check before you use
        numTE = NX(8); 
        id = 1:numTE:length(Lin);
        Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); id=[];
    end
    idxToUse = [ Lin(1) , Par(1) ];%Make sure we don't extract information from outside the elliptical shutter
    
    %CALC GEOM FOR HIGH-RES DATASET
    %SCALING
    rec.Par.Mine.Asca=diag([rec.Enc.x.AcqDelta 1]);

    %ROTATION
    rec.Par.Mine.Arot = double( eye(4));
    rec.Par.Mine.Arot(1:3,1:3)=double( cat(2, rec.ismrmrd.hdr.buffer.read_dir(:,idxToUse(1),idxToUse(2)), rec.ismrmrd.hdr.buffer.phase_dir(:,idxToUse(1),idxToUse(2)), rec.ismrmrd.hdr.buffer.slice_dir(:,idxToUse(1),idxToUse(2)))); %RO-PE-SL to PCS
    Arot_tmp = rec.Par.Mine.Arot;
    rec.Par.Mine.Arot(1:3,1) = Arot_tmp(1:3,2);rec.Par.Mine.Arot(1:3,2) = Arot_tmp(1:3,1); 
    clear Arot_tmp

     %TRANSLATION
    for i=1:hdrUPLong_l
        if strcmp(rec.ismrmrd.hdr.hdrUPLong(i).name,'lGlobalTablePosTra_zihan')
            lGlobalTablePosTra = rec.ismrmrd.hdr.hdrUPLong(i).value;
        end
    end
    rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = double(rec.ismrmrd.hdr.buffer.position(1:3,idxToUse(1),idxToUse(2)));%RO-PE-SL to PCS
    if exist('lGlobalTablePosTra')
        rec.Par.Mine.Atra(3,4) = rec.Par.Mine.Atra(3,4) + lGlobalTablePosTra; % considered the position of table, by ZN
    end

    % account for fact that tranRaw is referred to the center of FOV, not the first element in the array
    N = rec.Enc.x.AcqN;%Set inf to make sure this element is replaced
    orig = ( ceil((N(1:3)+1)/2) - [0 0 .5] )';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
    rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)...
                                - rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;

    %COMBINED MATRIX
    rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
    MT=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

    %Align the parameter 
    rec.Geom.x.patientPosition = rec.ismrmrd.hdr.hdrConnection.measurementInformation.patientPosition;
    rec.Geom.x.PCS2RAS=getPCS2RAS(rec.Geom.x.patientPosition);
    rec.Geom.x.MT=rec.Geom.x.PCS2RAS*MT;
    % ZN: for the next step, I actually doesn't know why we need translate
    [~,rec.Geom.x.MT] = mapNIIGeom([], rec.Geom.x.MT,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED

    % CALC THE GEOM FOR REFSCAN, ACCORDING TO THE GEOM OF HIGH-RES
    fprintf('Computing geometry for refscan/ACS.\n');
%     rec.Enc.S.AcqN = multDimSize(rec.ismrmrd.refdata,1:3);
    rec.Enc.S.AcqN = rec.Enc.S.AcqN([2 1 3]); % permute to [Lin RO Par]
    [rec.Enc.S.AcqDelta, rec.Geom.S.MT] = mapNIIGeom([], rec.Geom.x.MT,'resampling',[],rec.Enc.x.AcqN,rec.Enc.S.AcqN);%TO CHECK:AD HOC - REVERSE ENGINEERED

    % ZN: permute back to [RO Lin Par]
    % high-res
    rec.Enc.x.AcqN([1 2 3]) = rec.Enc.x.AcqN([2 1 3]);
    [rec.Geom.x.MS , rec.Geom.x.MT] = mapNIIGeom(rec.Enc.x.AcqDelta , rec.Geom.x.MT,'permute', [2 1 3]);
    rec.Enc.x.AcqDelta = rec.Geom.x.MS;
    % ACS/refscan
    rec.Enc.S.AcqN([1 2 3]) = rec.Enc.S.AcqN([2 1 3]);
    [rec.Geom.S.MS, rec.Geom.S.MT] = mapNIIGeom(rec.Enc.S.AcqDelta, rec.Geom.S.MT,'permute', [2 1 3]);
    rec.Enc.S.AcqDelta = rec.Geom.S.MS;

end
%%%%%%%%%%%%%%%%%%%%%%%INVERSION%%%%%%%%%%%%%%%%%%%%%%%
%INVERT NOISE
if 0 %isfield(TW,'noise') % ZN: disable noise-prewhitening for gadgetron pipeline 
    %READ NOISE
    if rec.Alg.removeOversamplingReadout;TW.noise.flagRemoveOS = true;end
    rec.N.(t)=TW.noise.unsorted;
    if gpu;rec.N.(t)=gpuArray(rec.N.(t));end    
    %INVERT READOUT
    rec.N.(t)=gather(invertReadout(rec.N.(t),rec.Enc.(t).GridMatrix));
else
    rec.N.(t)=[];
end

%INVERT PHASE
if typ==3;fieldPhase='phasecor';
elseif typ==2;fieldPhase='refscanPC';
else error('Unrecognized field: %s',field);
end
if 0 %isfield(TW,fieldPhase) && ~isfield(rec,'PS') %ZN: this part is skipped for ref & high-res
    rec.P.(t)=dynInd(TW.(fieldPhase),{':'},ND);   
    %rec.P.(t)=dynInd(TW.refscanPC,{':'},ND);
    rec.P.(t)=sum(rec.P.(t),3);%Sum over the phase encodes, assumes that typically one, but sometimes with zeros 
    rec.P.(t)=sum(rec.P.(t),6)./(sum(single(rec.P.(t)~=0),6)+eps);%Sum over the averages (changed recently to include the eps, necessary in case fMRI_2023_07_04)
    rec.P.(t)=mean(rec.P.(t),[5 9]);%Average over the repeats       
    %rec.P.(t)=prod(rec.P.(t),2).^(1/size(rec.P.(t),2));   
    rec.P.(t)=mean(rec.P.(t),2); 
    if gpu;rec.P.(t)=gpuArray(rec.P.(t));end
    rec.P.(t)=invertReadout(rec.P.(t),rec.Enc.(t).GridMatrix);
%    visGhost(dynInd(rec.P.(t),1,9))
    rec.P.(t)=gather(ridgeDetection(rec.P.(t),1,32));  
%    visGhost(dynInd(rec.P.(t),1,9))
%1
else
    rec.P.(t)=[];
end
% ZN: since nyquist correction & mask generation of ky and mb are all
% disabled for both high-res and refscan -> there's no need to keep these
% two parameters then
% rec.P.echoOrder=TW.(field).Seg(TW.(field).Sli==1 & TW.(field).Rep==1);
% rec.P.corrLine=TW.(field).Lin(TW.(field).Sli==1 & TW.(field).Rep==1);


%INVERT READOUT DATA
rec.(t)=[];
blkSz=8;%Depends on available GPU memory
for s=1:blkSz:rec.Enc.(t).CoilN;vS=s:min(s+blkSz-1,rec.Enc.(t).CoilN);%Inversion by channel groups
    %EXTRACT
    x=dynInd(rec.TWD.(t),vS,2);
    if gpu;x=gpuArray(x);end        
    %INVERT READOUT
    x=invertReadout(x,rec.Enc.(t).GridMatrix); 
    %ASSIGN
    rec.(t)=cat(4,rec.(t),gather(x));
end
rec.TWD.(t)=[];%Releasing memory

%NYQUIST CORRECTION
% ZN: this usually skipped for ref & high-res
if 0 %~isempty(rec.P.(t)) && strcmp(field,'image') && ~isfield(rec,'PS') 
    if gpu;rec.(t)=gpuArray(rec.(t));rec.P.(t)=gpuArray(rec.P.(t));end
    for n=1:size(rec.P.(t),11)
        in=rec.P.corrLine(rec.P.echoOrder==n);
        in=unique(in);%Recent addition to prevent situations where indexes are repeated, a bit risky
        rec.(t)=dynInd(rec.(t),in,2,dynInd(rec.(t),in,2).*exp(-1i*angle(dynInd(rec.P.(t),n,11))));
    end
    rec.(t)=gather(rec.(t));rec.P.(t)=gather(rec.P.(t));
end

%GENERATE MASKS FOR KY AND MB
rec.Ay=sum(abs(rec.(t)).^2,setdiff(1:ND,2));
rec.Ay=single(rec.Ay>1e-12);
if typ==3;rec.Ay=rec.Ay+(1-rec.Enc.(t).APF{2});end
if 0 %(~isempty(rec.P.(t)) || isfield(rec,'PS')) && strcmp(field,'image') % ZN: this is usually skipped for ref & high-res
    for n=1:max(rec.P.echoOrder)
        in=rec.P.corrLine(rec.P.echoOrder==n);   
        in=unique(in);%Recent addition to prevent situations where indexes are repeated, a bit risky
        rec.Ay(in)=rec.Ay(in)*n;
    end
end

%MULTIBAND
% ZN: I think this part skipped for both ref & high-res, but worth
% double-check
rec.Enc.(t).MultiBandFactor=1; % ZN: set as default
% if rec.Enc.(t).MultiBandFactor>1 && strcmp(field,'image')
%     MB=rec.Enc.(t).MultiBandFactor;
%     if MB==2;shift=2;elseif MB==3;shift=2;end
%     %bl=generateGrid(MB,[],MB,1);bl=bl{1}';
%     bl=generateGrid(shift,[],MB,1);bl=bl{1}';
%     NKYFull=length(rec.Ay);NKYSamp=sum(rec.Ay>0);
%     bl=bl(mod(0:NKYSamp-1,2)+1);%Blipping pattern
%     rec.Az=rec.Ay;rec.Az(:)=0;rec.Az(rec.Ay>0)=bl;
%     %sortedSliceOffsets=dynInd(rec.Enc.(t).SliceOffsets,rec.Enc.(t).SlicePositions,2,rec.Enc.(t).SliceOffsets);
%     NSl=length(rec.Enc.(t).SlicePositions);
%     sliceOffsets=resPop(repelem((0:MB-1)*2*pi/MB,NSl/MB),2,[],3);
%     if MB==2;slPhase=[0 0];elseif MB==3;slPhase=[0 0 pi];end;slPhase=resPop(repelem(slPhase,[NSl/MB]),2,[],3);
%     rec.Az=exp(1i*(rec.Az.*sliceOffsets+slPhase));
%     rec.Az=reshape(rec.Az,[1 NKYFull size(rec.Az,3)/MB MB]);    
%     rec.Az=ifftshiftGPU(rec.Az,2);
% else
%     rec.Az=[];
% end
rec.Az = [];
rec.Ay=ifftshiftGPU(rec.Ay,2);

%INVERT
blkSz=blkSz*2;
for s=1:blkSz:rec.Enc.(t).CoilN;vS=s:min(s+blkSz-1,rec.Enc.(t).CoilN);%Inversion by channel groups
    %EXTRACT
    x=dynInd(rec.(t),vS,4);
    if gpu;x=gpuArray(x);end        
    %INVERT
    for n=2:3;x=ifftInvert(x,n);end    
    %ASSIGN 
    rec.(t)=dynInd(rec.(t),vS,4,gather(x));
end

%REORDER, MAP SLICES TO THIRD DIMENSION, AND REMOVE UNUSED SLICES FOR MB
perm=1:ND;perm([3 5])=[5 3];
if 0 %~rec.Enc.(t).ThreeD
    rec.(t)=dynInd(rec.(t),rec.Enc.(t).SlicePositions,5,rec.(t));
    rec.(t)=permute(rec.(t),perm);
end

%%% WRITE METADATA TO JSON FILE (added by ZN)
writeJSONFlag=1;%Hard-coded to 1
if writeJSONFlag>0
    recJSON = rec;if nargout<1;rec=[];end
    recJSON.NY=size(recJSON.(t));
    %Remove all big arrays
    if isfield(recJSON,'S'),recJSON.S = [];end
    if isfield(recJSON,'TW'),recJSON.TW = [];end
    if isfield(recJSON,'TWD'),recJSON.TWD = [];end
    if isfield(recJSON,'N'),recJSON.N = [];end
    if isfield(recJSON,'P'),recJSON.P = [];end
    if isfield(recJSON,'Ay'),recJSON.Ay = [];end
    
    if isfield(recJSON,'M'),recJSON.M = [];end
    if isfield(recJSON,'B'),recJSON.B = [];end
    if isfield(recJSON,'x'),recJSON.x = [];end
    if isfield(recJSON,'xS'),recJSON.xS = [];end
    if isfield(recJSON,'W'),recJSON.W = [];end
    if isfield(recJSON,'ismrmrd'),recJSON.ismrmrd = [];end
    
    %Save JSON
    savejson('',recJSON,sprintf('%s.json',[rec.Nam.caseIn,'/JSON_',t]));%writeJSON has specific fields 
    fprintf('JSON file saved \n');
end


%%%%%%%%%%%%%%%%%%%%%%%AUXILIARY FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%
function x=invertReadout(x,gridMatrix)
    NF=size(x,1);%Input size    

    %PERMUTE COILS TO 4TH DIMENSION
    perm=1:ND;perm(2:4)=[3:4 2];
    x=permute(x,perm);

    %REGRIDDING MATRIX
    if ~isempty(gridMatrix);A=gridMatrix;else A=eye(NF,'like',x);end

    %INVERSION MATRIX
    F=build1DFTM(NF,0,isa(x,'gpuArray'));
    F=ifftshift(fftshift(F,2),1);

    %REMOVE OVERSAMPLING AND CROP TO DEFINED AREA
    F=dynInd(resampling(F,rec.Enc.(t).AcqN(1),3),vr,1);
    
    %APPLY
    x=aplGPU(F*A,x,1);
end

function x=ifftInvert(x,m)
    NF=size(x,m);%Input size
    %m

    %INVERSION MATRIX
    F=build1DFTM(NF,0,isa(x,'gpuArray'));
    F=ifftshift(fftshift(F,2),1);
    %F=circshift(F,ceil((rec.Enc.(t).OutN(m)-rec.Enc.(t).RecN(m))/2),1);
    %APPLY
    x=aplGPU(F,x,m);
    %m
    %any(isnan(x(:)))
end

end
