function rec=mapGeometry(fil,typ,rec,Par)

%MAPGEOMETRY   Maps the geometries of a set of datasets into the geometry of another dataset
%   Y=MAPGEOMETRY(Y,MT,REC)
%   * FIL is the filename of the dataset that needs to be mapped
%   * TYP are the data types to be mapped
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), 
%   the .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)), the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * PAR is some info about the files to be read
%   * REC is the modified reconstruction structure with the datasets mapped
%

gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);
parS=rec.Alg.parS;
if (strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPNeoBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiBra')) && ismember(rec.Par.Mine.Modal,9:10)
    parS.nDilate=repmat(parS.nDilate,[1 3]);
    parS.nDilate([1 2 3])=0;
else
    parS.nDilate=0;
end

suff=cell(1,length(typ));
assert(all(typ<=rec.Plan.NTypes),'Input types (%s) not supported, only %d types supported',sprintf(' %d',typ),rec.Plan.NTypes);
for n=1:length(typ);suff{n}=strcat(rec.Plan.TypeNames{typ(n)},rec.Plan.Suff);end
[y,~,MT]=readNII(fil,suff,gpu);

%MAP THE COORDINATES BACK
MTtarget=rec.Par.Mine.APhiAcq;
MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];
NPE=size(MTtarget,3);
MTT=repmat(MTT,[1 1 NPE]);
for s=1:NPE;MTtarget(:,:,s)=MTT(:,:,s)*MTtarget(:,:,s);end

for n=1:length(y);datTyp=rec.Plan.Types{typ(n)};  
    for s=1:NPE;MTT(:,:,s)=MT{n}\MTtarget(:,:,s);end    
    %TODO, IMPLEMENT USING SINC INTERPOLATION. IT SHOULD BE DONE IN TWO 
    %STEPS, FIRST AN ORTHOGONALIZATION OF THE MATRIX AND CONSEQUENT 
    %ROTATION AND THEN A SEPARABLE RESAMPLING OPERATION    
    if typ(n)==7;Nout=size(rec.y,4);else Nout=1;end

    %%%HERE WE HAD USED PAR, BUT IT GAVE PROBLEMS
    if ismember(typ(n),[7 27]) && ~isempty(Par) && ismember(Par.Mine.Modal,9:10) && rec.Alg.UseSBSensi==3;Neig=rec.Alg.parE.NCV;else Neig=1;end%Number of coils/Number of eigenmaps  
    Nsource=size(y{n});Nsource(end+1:3)=1;Nsource=Nsource(1:3);
    Ndestin=rec.Enc.FOVSize; 
    if isempty(rec.Par.Encoding.KzRange)                        
        assert(length(unique(rec.Par.Encoding.NrFids))==1,'Different stacks / mixes acquired with different NrFids%s',sprintf(' %d',rec.Par.Encoding.NrFids));               
        Ndestin(3)=size(rec.Par.Labels.Thicknesses,1)*rec.Enc.MB;    
    end
    NdestinReal=Ndestin;
    NdestinReal(1:3)=NdestinReal(1:3)./abs(rec.Alg.OverDec(1:3));
    
    shift=(NdestinReal-Ndestin)/2;
    rdGrid=generateGrid(Ndestin,gpuIn,Ndestin,ones(1,3),shift);
    if isfield(rec.Enc,'ExcGrid')
        rdGrid{3}=permute(rec.Enc.ExcGrid,[1 3 2]);        
        if gpuIn;rdGrid{3}=gpuArray(rdGrid{3});end
    end
    [dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rdGrid{1}(:),rdGrid{2}(:),rdGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;
    destinGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');dGrid{4}=[];
    if gpuIn;MTT=gpuArray(MTT);end
        
    x=single(zeros([Ndestin Nout NPE Neig]));    
    if gpu;x=gpuArray(x);end     
    x=complex(x);
    sdGrid=generateGrid(Nsource,gpuIn,Nsource,ones(1,3));
    [sGrid{1},sGrid{2},sGrid{3}]=ndgrid(sdGrid{1}(:),sdGrid{2}(:),sdGrid{3}(:));   
    warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId');
    for s=1:NPE    
        destinGridS=MTT(:,:,s)*destinGrid;
        for m=1:3
            dGrid{m}=reshape(destinGridS(m,:),Ndestin);
            dGrid{m}(dGrid{m}<min(sGrid{m}(:)))=min(sGrid{m}(:));
            dGrid{m}(dGrid{m}>max(sGrid{m}(:)))=max(sGrid{m}(:));
        end
        for t=1:Neig
            for m=1:Nout
                seff=s;
                z=dynInd(y{n},(t-1)*size(y{n},4)/Neig+mod(m-1+(seff-1)*Nout,size(y{n},4)/Neig)+1,4);
                if gpuIn;z=gpuArray(z);end
                z=interpn(sGrid{1},sGrid{2},sGrid{3},z,dGrid{1},dGrid{2},dGrid{3},'linear',0);
                if ~gpu;z=gather(z);end
                x=dynInd(x,[m s t],[4 5 6],z);
            end
        end
    end%Very quick, spline not supported in GPU    
    warning('on','all');
    
    %%CHECK THE GEOMETRY
    %if any(isnan(x(:)));fprintf('Acquired FOV not encompassed by the SENSE reference FOV\n');rec.Fail=1;return;end

    %First case is soft masking
    if typ(n)==12;rec.M=x;else rec.(datTyp)=x;end
    y{n}=[];
    if typ(n)==7 && rec.Alg.NoiseStand && rec.Par.Mine.Modal~=2 && isfield(rec,'N');rec.S=standardizeCoils(rec.S,rec.N);end%Coil standardization
    if ismember(typ(n),[8 12])
        voxsiz=rec.Enc.AcqVoxelSize;voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
        if typ(n)==8 && (~((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10) || ~rec.Alg.UseSBSensi || (isempty(rec.Par.Mine.AdHocArray) || ~((rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1))))
            rec.M=refineMask(rec.M,parS,voxsiz);
        end
        if isfield(rec,'S');rec.M(multDimSum(abs(rec.S),[4 6])<1e-9)=0;end%Outside the calibrated area
        if typ(n)==8;rec.M=single(rec.M>0.5);end
    end    
    if typ(n)~=12
        rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,typ(n));
        rec.Dyn.Typ2Wri(typ(n))=0;
    end
end
