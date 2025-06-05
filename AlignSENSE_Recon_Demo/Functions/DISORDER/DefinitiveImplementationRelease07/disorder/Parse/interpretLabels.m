function LabelsNew=interpretLabels(metaData,Labels)

%INTERPRETLABELS   Interprets the Labels of a .raw file
%   LABELSNEW=INTERPRETLABELS(METADATA,LABELS)
%   * METADATA is the metadata information as obtained from the feedHeader function
%   * {LABELS} is a set of labels from RF for comparison
%   * LABELSNEW are the interpreted labels
%

if nargin<2;Labels=[];end

ReconLabels=metaData(1:metaData(3));
OtherLabels=metaData(metaData(3)+1:end);
InterpretedLabels=zeros(1,metaData(3));

blockID=ReconLabels(1);
threePointOneField=ReconLabels(2);

InterpretedLabels([1 2 3])=1;

n=4;
b=0;
while 1    
    if metaData(n)==blockID
        InterpretedLabels(n:n+2)=1;
        b=metaData(n+1)+1;
        indB{b}=n+3:n+3+metaData(n+2)-1;
        block{b}=metaData(indB{b})';    
        n=n+3+metaData(n+2);
        if n>metaData(3)+3
            assert(n==metaData(3)+4,'Inconsistent block partition and recon label sizes: %d %d',n-4,metaData(3));
            fprintf('Completete reading: %d blocks\n',b);
            break
        end    
    else
        fprintf('Incomplete reading: %d of %d entries (%d blocks)\n',n-3,metaData(3),b);  
        break
    end
end
LabelsNew=[];

%%BLOCK 1-GENERAL INFORMATION
LabelsNew.NumberOfMixes=block{1}(1);
LabelsNew.NumberOfEchoes(1:block{1}(1),1)=block{1}(2);
LabelsNew.NumberOfEncodingDimensions(1:block{1}(1),1)=block{1}(3)/sum(LabelsNew.NumberOfEchoes);
InterpretedLabels(indB{1}(1:3))=1;
co=4;cd=co+sum(LabelsNew.NumberOfEncodingDimensions);
LabelsNew.VoxelSizes=block{1}(co:cd-1);
if LabelsNew.NumberOfEncodingDimensions(1)~=3
    LabelsNew.VoxelSizes(end+1:end+LabelsNew.NumberOfMixes)=block{1}(cd);
    InterpretedLabels(indB{1}(co:cd))=1;
end
LabelsNew.VoxelSizes=reshape(LabelsNew.VoxelSizes,block{1}(2),[]);%Problems with dimensions of array for two mixes only

co=12;cd=co+block{1}(1)-1;
LabelsNew.RepetitionTime=block{1}(co:cd);
InterpretedLabels(indB{1}(co:cd))=1;
co=14;cd=co+block{1}(1)-1;
LabelsNew.FlipAngle=block{1}(co:cd);
InterpretedLabels(indB{1}(co:cd))=1;
LabelsNew.ZReconLength=block{1}(16);
LabelsNew.TFEfactor=block{1}(20);
LabelsNew.EPIFactor=block{1}(35);
LabelsNew.DiffusionEchoTime=block{1}(38);
LabelsNew.DiffusionValues=block{1}(39);
LabelsNew.GradientOris=block{1}(40);
LabelsNew.FieldStrength=block{1}(42);
LabelsNew.ResonanceFreq=block{1}(43);
LabelsNew.ScanDuration=block{1}(52);
InterpretedLabels(indB{1}([16 20 35 38:40 42:43 52]))=1;
%InterpretedLabels(indB{1}([50 51]))=1;%PROBLEMATIC FOR VISUALIZING!
InterpretedLabelsBlock{1}=InterpretedLabels(indB{1});
%block{1}=block{1}.*(1-InterpretedLabelsBlock{1});
%figure;plot(block{1})%RESIDUAL PROBLEMS

%LabelsNew.ResonanceFreq
%Labels.ResonanceFreq

if ~isempty(Labels)%WE CHECK THE INFORMATION
    fiel=setdiff(fieldnames(LabelsNew),{'ResonanceFreq','VoxelSizes'});
    for n=1:length(fiel)      
        assert(all(size(LabelsNew.(fiel{n}))==size(Labels.(fiel{n}))),'Sizes not matched for field %s (new:%s / old:%s)',fiel{n},sprintf(' %d',size(LabelsNew.(fiel{n}))),sprintf(' %d',size(Labels.(fiel{n}))));
        assert(all(abs(LabelsNew.(fiel{n})-Labels.(fiel{n}))<1e-9),'Not matched element %s',fiel{n});
    end
end
    

%%%BLOCK 2-MIN K INFORMATION
LabelsNew.Kmin=block{2}(1:block{1}(3));
LabelsNew.Kmin=reshape(LabelsNew.Kmin,prod(block{1}(1:2)),[]);
LabelsNew.Kmin(:,end+1:3)=0;


%NONINTERPRETED FIELDS:
%                ASLNolabelTypes: 0
%                        ASLType: 'No'
%                AcquisitionMode: 'Cartesian'
%                     AdHocArray: [1x128 double]
%                      AngioMode: 'CE'
%                     Angulation: [0 0 0]
%               AverageOffcentre: [3.6455 -16.5552 0]
%                       CardSync: 'No'
%                        CoilNrs: []
%                  ConcomFactors: [1x384 double]
%                           Date: '20.07.2018'
%                      Diffusion: 0
%                    DiffusionAP: [1x128 double]
%               DiffusionBValues: [1x32 double]
%                    DiffusionFH: [1x128 double]
%                    DiffusionRL: [1x128 double]
%                    DynamicScan: 0
%                           Echo: 0
%                     FEARFactor: 0
%                            FOV: [375 148.8095 201.0000]
%                FastImagingMode: 'TFE'
%                    FatShiftDir: 'L'
%                       FlowComp: 0
%                    FoldOverDir: 'AP'
%                    GeoCorrPars: [1x354 single]
%             HeartPhaseInterval: 0
%                      KooshBall: 'no'
%                             Kt: 'No'
%                       KtFactor: 1
%                    KtReconMode: 'Blast'
%             KxOversampleFactor: 2
%                        KxRange: [-252 251]
%             KyOversampleFactor: 1.0300
%                        KyRange: [-20 50]
%             KzOversampleFactor: 1.2836
%                        KzRange: [-30 42]
%                  MPSOffcentres: [11 -2 0]
%                MPSOffcentresMM: [16.5552 -3.6455 0.9000]
%                            MTC: 0
%                     MinSliceNr: 0
%                            Mix: 0
%                      Multivenc: []
%                    NrInstances: 1
%                     NrSegments: 1
%                      NusEncNrs: [1x0 double]
%                      NusMethod: 0
%                     NusSamples: 0
%                      Offcentre: [3.6455 -16.5552 0]
%                    Orientation: 'TRA'
%                      PCAcqType: 'MPS'
%                     PDAFactors: [1x384 double]
%          PartialFourierFactors: [1 0.7000 0.8500]
%             PatientOrientation: 'Supine'
%                PatientPosition: 'HeadFirst'
%                      QuantFlow: 'No'
%                        RefScan: 'No'
%                           Rel4: 0
%                       RespComp: 'No'
%                       RespSync: 'Breathold'
%                    SENSEFactor: [1 2 1]
%                           SPIR: 1
%                     SampleSets: 1
%                        Samples: [252 99.0291 67]
%                   ScanDuration: 18.5992
%                       ScanMode: '3D'
%                  ScanTechnique: 'T1TFE'
%                       ScanType: 'Imaging'
%                      SliceGaps: 0
%                        Spectro: 0
%           SpiralLeadingSamples: 0
%                     StackIndex: 0
%                             TE: [1.4223 30 30]
%                    Thicknesses: [375 297.6190 258]
%                           Time: '18:49:29'
%                            UTE: 'no'
%                           Venc: [0 0 0]
%                            WFS: 0.6005
%                         XRange: [-128 127]
%                    XResolution: 256
%                         YRange: [-105 0]
%                    YResolution: 204
%                         ZRange: [-86 -1]
%                    ZResolution: 67
%                        Release: []
%                          Index: [1x1 struct]
%            OriginalLabelLength: 82946
%                CoilNrsPerStack: {[16x1 uint16]}