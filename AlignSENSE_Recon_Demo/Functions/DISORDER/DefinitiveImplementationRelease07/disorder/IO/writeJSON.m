function writeJSON(rec,dat,bids)

%WRITEJSON   Writes a json file with scan information
%   WRITEJSON(REC,DAT,{BIDS})
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   * DAT indicates the path and name of the file
%   * {BIDS} indicates whether to write under the bids format. It defaults
%   to 0. The BIDS fields are taken from https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/01-magnetic-resonance-imaging-data.html
%

if nargin<3 || isempty(bids);bids=0;end

Par.Encoding=rec.Par.Encoding;
Par.Labels=rec.Par.Labels;
Par.Scan=rec.Par.Scan;
Par.Mine=rec.Par.Mine;

if ~isfield(rec,'Enc');rec.Enc=[];end
if ~isempty(Par.Encoding)
    if isempty(Par.Encoding.KzRange)
        Par.Encoding.KzRange=[0 0];Par.Encoding.ZRange=[0 0];Par.Encoding.KzOversampling=0;
    end
    if isempty(Par.Encoding.ZRes)
        Par.Encoding.ZRes=0;Par.Encoding.ZReconRes=0;
    end

    if bids    
        if ~isfield(rec.Enc,'MB');rec.Enc.MB=1;end%We may want to look for this parameter, but not done for the moment
        codDir={'i','j','k'};
        N=size(rec.x);
        NewPar.Manufacturer='Philips Medical Systems';%Hardcoded from tag 0008,0070 of a DICOM
        NewPar.ManufacturersModelName='Achieva';%Tag 0008,1090 not present in our DICOMS
        %NewPar.DeviceSerialNumber='';%tag 0018,1000 is 38283 but anonymized following Anthony's suggestion
        %NewPar.StationName='';%tag 0008,1010 not present in DICOMS, we could set 'nnu' here but anonymized following Anthony's suggestion
        NewPar.SoftwareVersions='3.2.2\3.2.2.0';%Hardcoded from tag 0018,1020 in DICOM
        NewPar.MagneticFieldStrength=Par.Labels.FieldStrength/1000;
        NewPar.ReceiveCoilName='Dual coil';%Hardcoded from tag 0018,1250 in DICOM
        %NewPar.ReceiveCoilActiveElements='';???
        %NewPar.GradientSetType='';???
        %NewPar.MRTransmitCoilSequence='';%Should correspond to DICOM tag 0018,9049 but a dcmdump hasn't given anything meaningful for that field
        %%NewPar.MatrixCoilMode='';%Recommended only if used, it refers to coil combination (i.e. compression I think)
        NewPar.CoilCombinationMethod='SENSE';%We've opted for this denomination to point that the phase is used when combining the coils

        NewPar.PulseSequenceType=sprintf('%s %s %s %s',Par.Scan.Technique,Par.Scan.FastImgMode,Par.Scan.ScanMode,Par.Scan.AcqMode);
        %NewPar.ScanningSequence=;%In DICOM it would be 0018,0020 but we don't have direct acces to DICOMS
        %NewPar.SequenceVariant=;%In DICOM it would be 0018,0021 but we don't have direct acces to DICOMS
        %NewPar.ScanOptions=;%It would be 0018,0022 but not present in our DICOMS
        NewPar.SequenceName=Par.Scan.Technique;%It would be 0018,0024 but not present in our DICOMS, instead tag 0018,9005 is called PulseSequenceName which seems to correspond to Par.Scan.Technique
        NewPar.PulseSequenceDetails=sprintf('FastImgMode: %s / ScanMode: %s / AcqMode: %s',Par.Scan.FastImgMode,Par.Scan.ScanMode,Par.Scan.AcqMode);
        NewPar.NonlinearGradientCorrection=false;

        NewPar.NumberShots=0;
        if ismember(Par.Mine.Modal,2:3)
            NewPar.NumberShots=1;
        elseif Par.Labels.TFEfactor~=1 && ~strcmp(Par.Scan.ScanMode,'3D')
            NewPar.NumberShots=round((diff(Par.Encoding.KyRange(1,:))+1)/(Par.Labels.PartialFourierFactors(1,2)*Par.Labels.TFEfactor));
        elseif Par.Labels.TFEfactor~=1 && isfield(rec,'Assign')
            NewPar.NumberShots=round(length(rec.Assign.z{2})/Par.Labels.TFEfactor);
        elseif ~strcmp(Par.Scan.ScanMode,'3D')
            NewPar.NumberShots=round((diff(Par.Encoding.KyRange(1,:))+1)/(Par.Labels.PartialFourierFactors(1,2)*Par.Labels.ZReconLength));
        elseif isfield(rec,'Assign')
            NewPar.NumberShots=round(length(rec.Assign.z{2})/Par.Labels.ZReconLength);
        end
        NewPar.ParallelReductionFactorInPlane=prod(Par.Labels.SENSEFactor);%Corresponds to DICOM tag 0018,9069
        NewPar.ParallelAcquisitionTechnique='SENSE';%Corresponds to DICOM tag 0018,9078, not present in our DICOMS
        [NewPar.PartialFourier,pfd]=min(Par.Labels.PartialFourierFactors(1,:));%Corresponds to DICOM tag 0018,9081 which in our case is NO, but BIDS talks about fractions
        if NewPar.PartialFourier~=1;NewPar.PartialFourierDirection=codDir{pfd};end%Corresponds to DICOM tag 0018,9036 which in our case says FREQUENCY, but BIDS seems to use i,j,k for everything else
        if ~isempty(Par.Mine.diInfo)
            if Par.Mine.Nat==1;ind={'j','-j','i','-i'};
            elseif Par.Mine.Nat==2;ind={'-j','j','-i','i'};
            elseif Par.Mine.Nat==3;ind={'-i','i','j','-j'};
            else ind={'i','-i','-j','j'};
            end
            NewPar.PhaseEncodingDirection=[];
            if size(Par.Mine.diInfo,2)==6;diInfo=Par.Mine.diInfo(Par.Mine.diInfo(:,6)==1,:);else diInfo=Par.Mine.diInfo;end            
            for n=1:size(diInfo,1);NewPar.PhaseEncodingDirection=cat(2,NewPar.PhaseEncodingDirection,ind(diInfo(n,5)));end
            if NewPar.PartialFourier~=1;NewPar.PartialFourierDirection=NewPar.PhaseEncodingDirection;end
        elseif strcmp(Par.Scan.ScanMode,'3D')
            NewPar.PhaseEncodingDirection=codDir(2:3);
        else
            NewPar.PhaseEncodingDirection=codDir{2};
        end
        if ismember(Par.Mine.Modal,9:10) && isfield(rec.Enc,'kGrid')
            NPE=length(rec.Enc.kGrid{2});               
            if isfield(Par.Mine,'ES')
                ES=Par.Mine.ES;
            else    
                WFS=rec.Par.Labels.WFS;
                ES=WFS/(434.2902*(NPE+1));%This is described at reverseDistortion.m; there are three options: it is 427.8888 (standard Philips), 434.215 (Achieva scanners) / 434.2902 (PPE method)). Note it is still 427.8888 in reverseDistortion.m
            end
            ESFOV=ES*NPE/N(2);
            NewPar.EffectiveEchoSpacing=ESFOV;         
        else
            NewPar.EffectiveEchoSpacing=0;
        end
        NewPar.TotalReadoutTime=NewPar.EffectiveEchoSpacing*(N(2)-1);

        if Par.Mine.Modal~=3;NewPar.EchoTime=Par.Scan.TE/1000;else NewPar.EchoTime=Par.Labels.TE/1000;end%DICOM 0018,0081 in seconds
        NewPar.InversionTime=0;%DICOM 0018,0082 in seconds, not present in reconFrame, we set it to 0 and manually insert it for T1 sequences    
        if NewPar.NumberShots==1 && ~strcmp(Par.Scan.ScanMode,'3D')
            NewPar.SliceTiming=Par.Scan.TR(1)*rec.Enc.MB/(N(3)*1000);
            NewPar.SliceTiming=NewPar.SliceTiming*(0:N(3)/rec.Enc.MB-1);        
            if isfield(rec,'Assign') && numel(rec.Assign.z{8}(:))==N(3)/rec.Enc.MB;NewPar.SliceTiming(rec.Assign.z{8}(:)+1)=NewPar.SliceTiming;end
            NewPar.SliceTiming=repmat(NewPar.SliceTiming,[1 rec.Enc.MB]);
        end
        if ~strcmp(Par.Scan.ScanMode,'3D');NewPar.SliceEncodingDirection=codDir{3};end    
        if ~ismember(Par.Mine.Modal,9:10) && isfield(rec.Enc,'kGrid')
            NRO=length(rec.Enc.kGrid{1});    
            WFS=rec.Par.Labels.WFS;
            DT=WFS/(434.2902*(NRO+1));%This is described at reverseDistortion.m; there are three options: it is 427.8888 (standard Philips), 434.215 (Achieva scanners) / 434.2902 (PPE method)). Note it is still 427.8888 in reverseDistortion.m
            NewPar.DwellTime=DT;%Note it is not effective which would be obtained by dividing as in ESFOV
        end

        NewPar.FlipAngle=Par.Scan.FlipAngle;    
        %NewPar.NegativeContrast%Probably not very important
        if ~strcmp(Par.Scan.ScanMode,'3D');NewPar.MultibandAccelerationFactor=rec.Enc.MB;end

        %NewPar.AnatomicalLandmarkCoordinates%Probably not very important

        %NewPar.InstitutionName='King''s College London';%Tag 0008,0080 not present in our DICOMS, removed due to potential concerns with anonimity
        %NewPar.InstitutionAddress='St Thomas''s Hospital, Westminster Bridge Rd, Lambeth, London SE1 7EH';%Tag 0008,0081 not present in our DICOMS due to potential concerns with anonimity
        %NewPar.InstitutionalDepartmentName='Neonatal Dept.';%From DICOM tag 0008,1040, removed due to potential concerns with anonimity

        %if ismember(Par.Mine.Modal,5:7);NewPar.ContrastBolusIngredient;end%DICOM TAG 0018,1048, probably not very important

        NewPar.RepetitionTime=rec.Par.Scan.TR(1)/1000;%It was within fMRI but think it is important also in other modalities
        if ismember(Par.Mine.Modal,9)        
            %NewPar.VolumeTiming%Not necessary when using TR and delay time
            NewPar.TaskName='rest';

            NewPar.NumberOfVolumesDiscardedByScanner=0;
            NewPar.NumberOfVolumesDiscardedByUser=0;
            NewPar.DelayTime=0;
            %NewPar.AcquisitionDuration%Not necessary
            %NewPar.DelayAfterTrigger%Not necessary for resting state

            NewPar.Instructions='none';
            NewPar.TaskDescription='resting state sleep';
            NewPar.CogAtlasID='http://www.cognitiveatlas.org/id/trm_54e69c642d89b/';
            %NewPar.CogPOID%Unclear
        end

        NewPar.AcqVoxelSize=Par.Scan.AcqVoxelSize(1,:);%Added by me
        NewPar.RecVoxelSize=Par.Scan.AcqVoxelSize(1,:);%Added by me
        if ~strcmp(Par.Scan.ScanMode,'3D');NewPar.SliceGap=Par.Scan.SliceGap(1);end%Added by me

        if ismember(Par.Mine.Modal,10) && ~isempty(Par.Mine.diInfo)
            %PENDING SOME CHECKS WITH DAAN
            bvec=diInfo(:,1:3)';%FSL FORMAT, THREE ROWS AND ORIENTED AS THE IMAGE DATA
            %bvec([1 2,:)=-bvec([2 1],:);%It seems like this for neonates according to mrtrix conversion to fsl
            bvec(2,:)=-bvec(2,:);%It seems like this for neonates according to Matteo's email on 2019_02_05
            dlmwrite(sprintf('%s_Aq.bvec',dat),bvec,'delimiter',' ','precision','%.5f')
            bval=diInfo(:,4)';%FSL FORMAT, JUST CONCATENATING THEM
            dlmwrite(sprintf('%s_Aq.bval',dat),bval,'delimiter',' ','precision','%d');
        end
        savejson('',NewPar,sprintf('%s_AqBIDS.json',dat));
    else
        for s=1:2;NewPar.Encoding.KRange(s,:)=[Par.Encoding.KxRange(s) Par.Encoding.KyRange(s) Par.Encoding.KzRange(s)];end
        NewPar.Encoding.KOversampling=[Par.Encoding.KxOversampling(1) Par.Encoding.KyOversampling(1) Par.Encoding.KzOversampling(1)];
        NewPar.Encoding.Res=[Par.Encoding.XRes(1) Par.Encoding.YRes(1) Par.Encoding.ZRes(1)];
        NewPar.Encoding.ReconRes=[Par.Encoding.XReconRes(1) Par.Encoding.YReconRes(1) Par.Encoding.ZReconRes(1)];
        NewPar.Encoding.PartialFourierFactor=Par.Labels.PartialFourierFactors(1,:);
        NewPar.Encoding.SENSEFactor=Par.Labels.SENSEFactor;

        if isfield(Par.Labels,'AdHocArray')
            if Par.Labels.AdHocArray(1)==101
                NewPar.Encoding.MultibandFactor=Par.Labels.AdHocArray(4);
                NewPar.Encoding.MultibandShift=Par.Labels.AdHocArray(7);
                NewPar.Encoding.MultibandGap=Par.Labels.AdHocArray(2);
            else
                NewPar.Encoding.MultibandFactor=1;
                NewPar.Encoding.MultibandShift=0;
                NewPar.Encoding.MultibandGap=0;
            end
        end

        NewPar.Labels.AcquisitionMode=Par.Labels.AcquisitionMode;
        NewPar.Labels.EPIFactor=Par.Labels.EPIFactor;
        NewPar.Labels.FatShiftDir=Par.Labels.FatShiftDir;
        NewPar.Labels.FoldOverDir=Par.Labels.FoldOverDir;
        NewPar.Labels.FastImagingMode=Par.Labels.FastImagingMode;
        NewPar.Labels.FieldStrength=Par.Labels.FieldStrength;
        NewPar.Labels.FlipAngle=Par.Labels.FlipAngle;
        %NewPar.Labels.GeoCorrPars=Par.Labels.GeoCorrPars;
        NewPar.Labels.NumberOfEncodingDimensions=Par.Labels.NumberOfEncodingDimensions;
        %NewPar.Labels.Offcentre=Par.Labels.Offcentre;
        NewPar.Labels.Orientation=Par.Labels.Orientation;
        NewPar.Labels.PatientOrientation=Par.Labels.PatientOrientation;
        NewPar.Labels.PatientPosition=Par.Labels.PatientPosition;
        NewPar.Labels.RepetitionTime=Par.Labels.RepetitionTime;
        NewPar.Labels.Samples=Par.Labels.Samples;
        NewPar.Labels.ScanDuration=Par.Labels.ScanDuration;
        NewPar.Labels.ScanMode=Par.Labels.ScanMode;
        NewPar.Labels.ScanTechnique=Par.Labels.ScanTechnique;
        NewPar.Labels.EchoTime=Par.Labels.TE;
        NewPar.Labels.TFEfactor=Par.Labels.TFEfactor;
        NewPar.Labels.WFS=Par.Labels.WFS;
        NewPar.Labels.ReconstructionRelease=6;

        NewPar.Scan.AcqVoxelSize=Par.Scan.AcqVoxelSize;
        NewPar.Scan.RecVoxelSize=Par.Scan.RecVoxelSize;
        NewPar.Scan.SliceGap=Par.Scan.SliceGap(1);
        NewPar.Scan.FOV=Par.Labels.FOV;

        savejson('',NewPar,sprintf('%s_Aq.json',dat));
    end
end