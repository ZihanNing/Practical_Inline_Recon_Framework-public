function rec=reconstructData(rec)

%RECONSTRUCTDATA   Reconstructs data from acquired FOV to reconstructed FOV
%and estimates the parameters of the system matrix considered by the
%reconstruction model
%   REC=RECONSTRUCTDATA(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with reconstructed data and
%   surrogate information (rec.(rec.Plan.Types)) and removed information 
%   for correction (rec.Corr.(rec.Plan.Types)) and sorting 
%   (rec.Assign(rec.Plan.Types))
%

if rec.Fail;return;end
gpu=rec.Dyn.GPU;

name=rec.Names.Name;

%TO SAVE MEMORY
if rec.Dyn.GPU==2 && isfield(rec,'Recursive') && rec.Recursive==1
    for n=[7:8 27];datTyp=rec.Plan.Types{n};
        if isfield(rec,datTyp);rec.(datTyp)=gpuArray(rec.(datTyp));end
    end
end

if isfield(rec,'adding') && rec.adding==1;rec.Fail=3;return;end

%ESTIMATE/READ SENSITIVITY MAPS / B0 MAPS
if ((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10) && rec.Alg.DistoSensi;tyC={'S','B'};else tyC={'S'};end
if rec.Alg.UseSBSensi && ((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10)
    mo={[2 rec.Par.Mine.Modal],3};
    fo={{'Re-Se',numbe2Modal(rec.Par.Mine.Modal)},{'Re-F0'}};
    ex={'_Se','_B0'};
else
    mo={2,3};
    fo={{'Re-Se'},{'Re-F0'}};
    ex={'_Se','_B0'};
end
ty={7:8,11};
if rec.Alg.UseSoftMasking>0;ty{1}=cat(2,ty{1},12);end
if rec.Alg.UseSBSensi==3 && ~isempty(rec.Par.Mine.AdHocArray) && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1;ty{1}=cat(2,ty{1},27);end
for n=1:length(tyC)
    if ~isfield(rec,tyC{n})%Assign sensitivities        
        indREF=find(ismember(rec.Names.prot.A_Modal(1:rec.Names.ind-single(~(rec.Par.Mine.Modal==2))),mo{n}));        
        if ~isempty(indREF)%Candidate references
            indPrev=length(indREF);
            while 1%Check that is a valid reference          
                targetFile=strtrim(rec.Names.prot.B_FileName(indREF(indPrev),:));
                matFile=fullfile(rec.Names.headFolder,sprintf('%s.mat',targetFile));                 
                load(matFile);
                %THIS LINE DOES THE CHECK
                if ~(Par.Mine.Modal==9 && strcmp(Par.Scan.Technique,'SEEPI')) && ...
                    (length(Par.Parameter2Read.loca)>1 || length(Par.Parameter2Read.kz)>1) && ...
                    (isempty(Par.Mine.AdHocArray) || ~((Par.Mine.AdHocArray(1)==101 && Par.Mine.AdHocArray(4)~=1))) && ...
                    ((~isempty(rec.Par.Mine.AdHocArray) && rec.Par.Mine.AdHocArray(1)==101 && rec.Par.Mine.AdHocArray(4)~=1) || ismember(Par.Mine.Modal,[2 3])) && ...%Only FFE ref if it is SB 
                    (ismember(Par.Mine.Modal,[2 3]) || length(rec.Par.Mine.pedsUn)==length(Par.Mine.pedsUn) || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPEpiPha') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduPha'));break;end
                    %(ismember(Par.Mine.Modal,[2 3]) || length(rec.Par.Mine.pedsUn)==length(Par.Mine.pedsUn) || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodBra') || strcmp(rec.Par.Mine.curStud.Stud,'dHCPTodPha') || (strcmp(rec.Par.Mine.curStud.Stud,'dHCPAduBra') && strcmp(rec.Par.Mine.curStud.IssuId,'DiscardDWI')));break;end
                indPrev=indPrev-1;
                assert(indPrev~=0,'No valid reference %s for dataset %s',numbe2Modal(mo{n}(1)),name);
            end
            indFo=find(mo{n}==rec.Names.prot.A_Modal(indREF(indPrev)),1);
            fileSENSE=fullfile(rec.Names.pathOu,fo{n}{indFo},targetFile,'');            
            if indREF(indPrev)~=rec.Names.ind && (rec.Alg.ReestSensi || ~exist(strcat(fileSENSE,ex{n},rec.Plan.Suff,'.nii'),'file'))%We call the reference reconstruction
                recS.Names.Name=targetFile;
                matFile=fullfile(rec.Names.headFolder,sprintf('%s.mat',targetFile));rawFile=fullfile(rec.Names.pathIn,sprintf('%s.raw',targetFile));
                recS.Names.rawFile=rawFile;recS.Names.matFile=matFile;recS.Names.pathOu=rec.Names.pathOu;recS.Names.pathIn=rec.Names.pathIn;recS.Names.prot=rec.Names.prot;recS.Names.ind=indREF(indPrev);recS.Names.headFolder=rec.Names.headFolder;
                assert(logical(exist(rawFile,'file')),'Could not find the .raw file %s for sensitivity map estimation for file %s',rawFile,name);            
                load(matFile);
                recS.Par=Par;Par=[];
                %recS.Recursive=1;
                recS=reconPipeline(recS);               
                if recS.Fail;fprintf('Problem reconstructing the reference. See previous message\n');rec.Fail=1;return;end  
                recS=[];
            end
            if exist(strcat(fileSENSE,ex{n},rec.Plan.Suff,'.nii'),'file') && indREF(end)~=rec.Names.ind%We read the reference reconstruction
                assert(logical(exist(strcat(fileSENSE,ex{n},rec.Plan.Suff,'.nii'),'file')),'Sensitivity .nii file not found');     
                rec=mapGeometry(fileSENSE,ty{n},rec,Par);
            elseif mo{n}(indFo)==2%We solve the reference reconstruction
                rec=solveS(rec);%We perform the sensitivity estimation
            end
        end
    end
end

if isfield(rec,'B') && rec.Dyn.FirstIt
    rec.u=rec.S;rec.B=-rec.B;%Before I was negating rec.B but I think is not necessary, for fMRI seems to work better, for DWI seems to work worse... 
    rec.Dyn.Typ2Wri(11)=1;
    rec=reverseDistortion(rec,2.3,1,0);
    rec.Dyn.Typ2Wri(11)=0;   
    rec.S=rec.u;rec.B=-rec.B;
    rec=rmfield(rec,'u');
end

%REPEATED SCANS FOR MODALITY 7
if isfield(rec,'addRec')
    rec.y=gather(rec.y);
    for n=1:length(rec.addRec)
        rec.addRec{n}.adding=1;
        rec.addRec{n}=reconPipeline(rec.addRec{n});
        if rec.addRec{n}.Fail==1;rec.Fail=1;return;end
        rec.y=cat(5,rec.y,gather(rec.addRec{n}.y));
        for d=2:3;rec.Assign.z{d}=cat(d,rec.Assign.z{d},rec.addRec{n}.Assign.z{d});end       
        if isfield(rec.addRec{n},'y');rec.addRec{n}=rmfield(rec.addRec{n},'y');end
    end     
end

%RECONSTRUCT
if rec.Par.Mine.Modal~=2    
    if ismember(rec.Alg.SaveRaw,1:2);writeRaw(rec);end
    
    if ~isempty(rec.Alg.RemoveChannels)
        rch=rec.Alg.RemoveChannels;mch=setdiff(1:size(rec.y,4),rch);
        yr=dynInd(rec.y,rch,4);Sr=dynInd(rec.S,rch,4);
        rec.y=dynInd(rec.y,mch,4);rec.S=dynInd(rec.S,mch,4);
    end
    
    if rec.Alg.AlignedRec
        if rec.Par.Mine.Modal==7
            %PREVENT MEMORY OVERLOAD
            NX=size(rec.M);
            if prod(NX(1:3))>rec.Dyn.MaxMem(2)
                typ2Rec=rec.Dyn.Typ2Rec;
                for n=typ2Rec'        
                    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};rec.(datTyp)=gather(rec.(datTyp));end
                end
            end            
            
            %nam='GT';%nam='Q1';%nam='Q2';%nam='Q3';%nam='MPRAGE';%nam='SPGR';%nam='FLAIR';%nam='TSE';%nam='BSSFP';
            %save(sprintf('/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/DISORDERPaper/MatlabStructures/rec%s.mat',nam),'rec','-v7.3');
            %1
            %pause
            rec=solveXT(rec);
        else
            rec=solveXTMSOld(rec);
        end
    else        
        rec=solveX(rec);
    end
    
    if ~isempty(rec.Alg.RemoveChannels)
        rec.y=cat(4,rec.y,yr);rec.S=cat(4,rec.S,Sr);
        rec.y=dynInd(rec.y,[mch rch],4,rec.y);rec.S=dynInd(rec.S,[mch rch],4,rec.S);
    end
    
    if ismember(rec.Alg.SaveRaw,3:4);writeRaw(rec);end
end
if rec.Fail;return;end

%ESTIMATE B0 FIELD
if rec.Par.Mine.Modal==3;rec=solveB0(rec);end

%SOLVE FOR SENSITIVITIES FOR EPI-MB UNFOLDING
if ((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10) && rec.Alg.UseSBSensi>0 && rec.Enc.RecSize(3)==rec.Enc.FOVSize(3) && ~isfield(rec,'Recursive')
    if ~isfield(rec,'yProv');rec.yProv=rec.y;elseif rec.Par.Mine.Modal==10;rec.yProv=cat(11,rec.yProv,rec.y);else rec.yProv=cat(5,rec.yProv,rec.y);end
    if rec.Dyn.LastIt
        %CONCATENATE AND REVERSE MARGOSIAN
        rec.y=rec.yProv;rec=rmfield(rec,'yProv');
        if rec.Alg.MargosianFilter%Undoing, important for masking later
            Enc=rec.Enc;
            overDec=rec.Alg.OverDec;
            overDec(overDec<0)=1;
            Enc.FOVSize=round(Enc.FOVSize./overDec);
            rec.x=margosianFilter(rec.x,Enc,1,0);
        end
        
        %CALL THE SENSITIVITY ESTIMATION
        M=rec.M;%RECENT CHANGE
        %save('/home/lcg13/Work/DataDefinitiveImplementationDebug07/rec.mat','rec','-v7.3');
        %1
        %pause
        if rec.Alg.UseSBSensi==2;rec=solveS(rec);elseif ismember(rec.Alg.UseSBSensi,3:4);rec=solveESPIRIT(rec);end
        rec.M=M;M=[];%RECENT CHANGE
        
        %CALL THE RECONSTRUCTION AGAIN
        recSB.Names=rec.Names;
        assert(logical(exist(recSB.Names.rawFile,'file')),'Could not find the .raw file %s for SB reconstruction for file %s',recSB.Names.rawFile,name);            
        load(recSB.Names.matFile,'Par');
        recSB.Par=Par;Par=[];
        recSB.Recursive=1;
        tyN=[7 8];
        if rec.Alg.UseSBSensi==3;tyN=cat(2,tyN,27);end
        for n=tyN;datTyp=rec.Plan.Types{n};recSB.(datTyp)=rec.(datTyp);end

        for n=[3 5 6 7 8 9 12 27];datTyp=rec.Plan.Types{n};
            if isfield(rec,datTyp);rec.(datTyp)=gather(rec.(datTyp));end 
        end  
        if rec.Dyn.GPU==2            
            for n=tyN;datTyp=rec.Plan.Types{n};
                recSB.(datTyp)=gather(recSB.(datTyp));
            end
        end
        recSB=reconPipeline(recSB);
        if recSB.Fail;fprintf('Problem reconstructing the SB reference. See previous message\n');rec.Fail=1;return;end  
        recSB=[];
        if gpu
            for n=[3 5 6 7 8 9 12 27];datTyp=rec.Plan.Types{n};
                if isfield(rec,datTyp);rec.(datTyp)=gpuArray(rec.(datTyp));end
            end
        end
        rec.Dyn.Typ2Wri(:)=0;
    end
end
if rec.Dyn.LastIt && isfield(rec,'Recursive') && ((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10) && rec.Alg.UseSBSensi>0  
    tyN=[7 8];
    if rec.Alg.UseSBSensi==3;tyN=cat(2,tyN,27);end
    for n=tyN
        if ~any(ismember(rec.Dyn.Typ2Rec,n));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end
    end
end

rec=rmfield(rec,'y');rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==5 | rec.Dyn.Typ2Rec==6)=[];
if ~any(ismember(rec.Dyn.Typ2Rec,8));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,8);end
%if ~any(ismember(rec.Dyn.Typ2Rec,27));rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,27);end
if ~any(ismember(rec.Dyn.Typ2Rec,9)) && rec.Alg.EstimGFact(1) && ~rec.Alg.AlignedRec;rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,9);end
if ~any(ismember(rec.Dyn.Typ2Rec,7)) && rec.Alg.WriteSensi;rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,7);end
if rec.Alg.WriteSensi;rec.Dyn.Typ2Wri(7:8)=1;end

for n=[7 8 9 10 12 16 27];datTyp=rec.Plan.Types{n};
    if rec.Dyn.Debug>=2;printResults;end
end

function printResults
    if isfield(rec,datTyp) && ~isempty(rec.(datTyp));fprintf('Size of data typ %s:%s\n',datTyp,sprintf(' %d',size(rec.(datTyp))));end
end

end