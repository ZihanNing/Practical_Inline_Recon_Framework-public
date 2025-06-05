function rec=reconPlanner(rec)

%RECONPLANNER   Plans the memory management and desired output of a given 
%reconstruction and generates new structures for management
%   REC=RECONPLANNER(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (.Names), the status of the reconstruction (.Fail) 
%   and the .lab information (.Par)
%   ** REC is a reconstruction structure with filled plan information 
%   (.Plan for fixed plan information and .Dyn for dynamic plan 
%   information)
%

Plan.Quick=0;%For quick (1) and ultra-quick (2) reconstructions, mainly for testing DWI

name=rec.Names.Name;

%INDICATES SUFFIX
%REMOVE NOPF
Plan.Suff='';%Suffix to be applied when reading and writing data
Plan.SuffOu='';%Suffix to be applied only when writing data

%INDICATES WHETHER TO PRINT DEBUG INFORMATION(0->3 IN LEVELS OF VERBOSITY NOTHING/TIMES/BASIC/ALL)
%Dyn.Debug=2;
Dyn.Debug=0;

Plan.Dims={'size','ky','kz','chan','dyn','card','echo','loca','mix','extr1','extr2','aver'};Plan.NDims=length(Plan.Dims);
%1-> Raw / 2-> Empty / 3-> EPIRead / 4-> Navigator / 5-> Noise / 6-> Spectra
%7-> Sensitivities / 8-> Masks / 9-> G-facts / 10-> Chi2-maps / 11-> B0 / 12-> Reconstruction
%13-> Undistorted / 14-> Per-volume alignment / 15-> Per-excitation alignment / 16-> Volumetric alignment / 17-> Motion transform / 18-> Filtered reconstruction
%19-> Noise level / 20-> Number of components / 21-> Asymptotic error / 22-> Residuals / 23-> Frequency stabilization data / 24-> Image after frequency stabilization
%25-> Tracking information / 26->Tracked data / 27 -> SensitivityEigenMaps / 28 -> Generic container for mapping to material coordinates / 29 -> Linear B0 varying with rotation / 30 -> spin history information
Plan.Types={'z','','P','C','N','y', ...
            'S','M','G','X','B','x', ...
            'u','v','e','d','T','r', ...
            's','p','a','n','F','w', ...
            't','b','W','E','D','H'};
Plan.TypeNames={'Ra','','Ny','Na','Ga','Sp', ...
                'Se','Ma','No','Ch','B0','Aq', ...
                'Un','Vo','Ex','Di','Tr','Re', ...
                'Si','Co','Ae','Er','Fr','Ws', ...
                'Fo','Vt','Ei','Ec','Bd','Sh'};Plan.NTypes=length(Plan.TypeNames);
Plan.CorrTypes={'random_phase','pda_fac','pda_index','meas_phase','sign'};Plan.NCorrs=length(Plan.CorrTypes);

%DECIDE ON WHAT VOLUMES ARE TO BE RECONSTRUCTED AND START GENERATING THE 
%DATA SIZES FROM THEM
Plan.Par2Rec=cell(1,Plan.NDims);%Echo/Dyn/Extr2
if isfield(rec,'recV');rec.Par.Parameter2Read.(Plan.Dims{11})=(0:find(rec.Par.Mine.splitInfo~=0,1,'last')-1)';end%For split scans
for n=2:Plan.NDims;Plan.Par2Rec{n}=single(rec.Par.Parameter2Read.(Plan.Dims{n}));end
Plan.Typ2Rec=single(rec.Par.Parameter2Read.typ);
Plan.DimsOutLoop=[5 6 9 10 11 12];%These dimensions are assumed to be always outside the loop
Plan.DimsInsLoop=[2 3 7 8];%These dimensions are assumed to be always inside the loop
%The user may uncomment these to perform quick reconstructions of a few volumes
%Plan.Typ2Rec=[1;5];
%Plan.Par2Rec{8}=9;
%Plan.Par2Rec{10}=0;
if Plan.Quick>=2;Plan.Par2Rec{5}=(0:19)';end%(0:19)';%(0:99)';%Dyn
%Plan.Par2Rec{6}=(0)';%Card
%Plan.Par2Rec{7}=(0)';%Echo
%Plan.Par2Rec{9}=1;%Mix
if Plan.Quick>=2;Plan.Par2Rec{11}=(0:9)';end%(0:9)';%(0:24)';%b-val
%if Plan.Quick>=2;Plan.Par2Rec{11}=(10)';end%(0:9)';%(0:24)';%b-val
%Plan.Par2Rec{11}=(20:24)';

for n=5:12;Plan.Par2Rec{n}(Plan.Par2Rec{n}>max(single(rec.Par.Parameter2Read.(Plan.Dims{n}))))=[];end
Plan.ReturnHeader=0;%To return the header information

for n=2:Plan.NDims
    %NOW WE ARE NOT CHECKING WHETER PAR2REC IS NOT CONTIGUOUS FOR THE COILS
    if n~=4 && ((n~=3 || isempty(rec.Par.Mine.AdHocArray) || rec.Par.Mine.AdHocArray(1)~=101) && any(diff(single(Plan.Par2Rec{n}))~=1));fprintf('Currently the Par2Rec %d has to be contiguous and it is%s for file %s\n',n,sprintf(' %d',Plan.Par2Rec{n}),name);rec.Plan=Plan;rec.Dyn=Dyn;rec.Fail=1;return;end
end
for n=[5 6 7 9 11];Dyn.Size(n)=length(Plan.Par2Rec{n});end%These may admit changes

%DECIDE ON THE BLOCK SIZES TO SPLIT THE RECONSTRUCTION
Dyn.Block2Rec=inf(1,Plan.NDims);%Echo/Dyn/Extr2
Dyn.Block2Rec(5)=10;%Dyn
Dyn.Block2Rec(6)=1e6;%Card
Dyn.Block2Rec(7)=1e6;%Echo
Dyn.Block2Rec(9)=1;%Mix-This is strictly required to push the mixes outside the loop
Dyn.Block2Rec(11)=5;%b-val

%INDEXES TO BE RECONSTRUCTED/READ AT A GIVEN ITERATION
Dyn.Ind2Rec=cell(1,Plan.NDims);
for n=2:Plan.NDims;Dyn.Ind2Rec{n}=Plan.Par2Rec{n};end
Dyn.Typ2Rec=rec.Par.Parameter2Read.typ;Dyn.Typ2Wri=ones(1,Plan.NTypes);Dyn.Typ2Wri([1:6 17 23 25])=0;

%PREFERABLE TO INCLUDE THE INDEXES TO THE ARRAY / CORRECTIONS AS A NEW ASSIGNMENT / Correction STRUCTURE
Plan.Assign=cell(1,Plan.NDims);for n=2:Plan.NDims;Plan.Assign{n}=single(rec.Par.Labels.Index.(Plan.Dims{n}));rec.Par.Labels.Index=rmfield(rec.Par.Labels.Index,Plan.Dims{n});end
Plan.Corr=cell(1,Plan.NCorrs);for n=1:Plan.NCorrs;Plan.Corr{n}=rec.Par.Labels.Index.(Plan.CorrTypes{n});rec.Par.Labels.Index=rmfield(rec.Par.Labels.Index,Plan.CorrTypes{n});end
Plan.Typ=rec.Par.Labels.Index.typ;rec.Par.Labels.Index=rmfield(rec.Par.Labels.Index,'typ');
Plan.MinK=cell(1,3);
for n=2:3;Plan.MinK{n}=min(Plan.Assign{n});end

%INDICATES WHETHER WE ARE IN THE FIRST/LAST INNER ITERATION OF THE 
%RECONSTRUCTION CHUNKS IN GENERAL AND FOR EACH DIMENSION AND THE BATCH WE
%ARE PROCESSING
Dyn.FirstIt=1;Dyn.FirstItDim=ones(1,Plan.NDims);
Dyn.LastIt=0;Dyn.LastItDim=ones(1,Plan.NDims);Dyn.LastItDim([5 6 7 9 11])=0;
Dyn.CurItDim=zeros(1,Plan.NDims);
Dyn.Batch=0;

%INDICATES WHETHER TO USE GPU PROCESSING
Dyn.BlockGPU=blockGPU;%To completely block GPU usage
Dyn.GPU=~Dyn.BlockGPU && gpuDeviceCount;
Dyn.MaxMem=[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
Siz(1)=max(abs(diff(rec.Par.Encoding.KxRange,1,2))+1,[],1)/min(rec.Par.Encoding.KxOversampling);
Siz(2)=max(abs(diff(dynInd(rec.Par.Encoding.YRange,1,3),1,2))+1,[],1);
Siz(3)=max(length(Plan.Par2Rec{3}),length(Plan.Par2Rec{8}));
Siz(4)=max(max(min(length(Plan.Par2Rec{5}),Dyn.Block2Rec(5)),min(length(Plan.Par2Rec{11}),Dyn.Block2Rec(11))),min(length(Plan.Par2Rec{12}),Dyn.Block2Rec(12)));
if rec.Par.Mine.Modal==2;Siz(4)=1;end

if ~isempty(rec.Par.Encoding.ZRange);Siz(3)=max(Siz(3),max(abs(diff(rec.Par.Encoding.ZRange,1,2))+1,[],1));end
if prod(Siz)>Dyn.MaxMem(1);Dyn.GPU=0;end
if prod(Siz)>Dyn.MaxMem(3) && Dyn.GPU;Dyn.GPU=2;end
Dyn.ExpSiz=Siz;%Expected size

%CREATE CONTAINER FOR RECONSTRUCTION OBJECT AND ASSIGN PLAN AND DYN FIELDS
rec.x=[];%NOTE THIS MAY BE NOT NECESSARY
rec.Plan=Plan;rec.Dyn=Dyn;
