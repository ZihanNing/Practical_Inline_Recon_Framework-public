function recOu=readData(rec)

%READDATA   Reads all or a specific piece of raw data
%   REC=READDATA(REC,{RAWFILE})
%   * REC is a reconstruction structure. At this stage it may contain the naming information (.Names), the status of the reconstruction (.Fail), the .lab information (.Par), 
%   the fixed plan information (.Plan) and the dynamic plan information (.Dyn) 
%   * RECOU is a reconstruction structure cell with filled data information
%   ((.Plan.Types)) for each of the channel types
%

if rec.Fail;recOu{1}=rec;return;end

%PROBLEMS READING DATA IN INGENIA-NNU DUE PROBABLY TO DSTREAM TECHNOLOGY AND COIL ARRRANGEMENT?

%DO SOME PREASSIGNMENTS
if ~isfield(rec.Names,'rawFile') || isempty(rec.Names.rawFile);rec.Names.rawFile=rec.Par.Filename.Data;end
rawFile=rec.Names.rawFile;

offset=rec.Par.Labels.Index.offset;
siz=rec.Par.Labels.Index.size;
name=rec.Names.Name;typ=rec.Plan.Typ;%For convenience

NC=length(offset)/length(typ);%In case the indexes are compressed through the different channels
chancompi=(1:NC)';
notSup=[2 4];
if rec.Par.Mine.Modal==10;notSup=2;end%Not supported data types
for n=notSup
    if any(rec.Plan.Typ2Rec==n);fprintf('Raw data contains typ %d which is currently not supported by the reconstruction\n',n);end
    rec.Plan.Typ2Rec(rec.Plan.Typ2Rec==n)=[];
end
rec.Dyn.Typ2Rec=rec.Plan.Typ2Rec;
if ~rec.Dyn.FirstIt;rec.Dyn.Typ2Rec(ismember(rec.Dyn.Typ2Rec,4:5))=[];end%Noise/Navigation is already read
if ~rec.Dyn.CurItDim(7);rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==3)=[];end%Ghosting calibration is already read
typ2Rec=rec.Dyn.Typ2Rec;


%READ THE DATA
if isempty(strfind(rec.Names.pathIn,'raw-ingenia'))%AS WITH INGENIA NO LUCK IN READING LIKE THIS
    [fidraw,su]= fopen(rawFile,'r','ieee-le');assert(fidraw>0,'Error %s for file %s',su,rawFile);
    for m=typ2Rec';datTyp=rec.Plan.Types{m};
        if isfield(rec,datTyp);rec=rmfield(rec,datTyp);end;rec.(datTyp)=[];
        indRead=(rec.Plan.Typ==m);
        if m~=5%We always read all the noise samples if indicated by typ2Rec
            for n=2:rec.Plan.NDims            
                if n~=4 || NC==1;indRead=(indRead & ismember(rec.Plan.Assign{n},rec.Dyn.Ind2Rec{n}));end
            end
        end
        indRead=find(indRead);        
        if ~isempty(indRead)
            blockLimits=find(diff(indRead)~=1);
            indReadFull=bsxfun(@plus,NC*(indRead-1)',chancompi);
            rec.Corr.(datTyp)=cell(1,rec.Plan.NCorrs);
            for n=1:rec.Plan.NCorrs
                if n~=2;rec.Corr.(datTyp){n}=rec.Plan.Corr{n}(indRead);else rec.Corr.(datTyp){n}=rec.Plan.Corr{n}(indReadFull(:));end
            end
            rec.Assign.(datTyp)=cell(1,rec.Plan.NDims);
            for n=2:rec.Plan.NDims
                if n~=4 || NC==1;rec.Assign.(datTyp){n}=rec.Plan.Assign{n}(indRead);end
            end            
            if ~isempty(blockLimits) && rec.Dyn.Debug>=1;fprintf('The data of typ %d forms %d blocks for file %s\n',m,length(blockLimits)+1,name);end
            blocks=[0;blockLimits;length(indRead)];
            for n=1:length(blocks)-1
                indReadBlock=indRead(blocks(n)+1:blocks(n+1));
                indReadBlock=bsxfun(@plus,NC*(indReadBlock-1)',chancompi);          
                offsetRead=offset(indReadBlock(:));
                sizRead=siz(indReadBlock(:));
                %In the future we should implement the functionality for
                %non-linearly growing offsets, such as for non consecutive
                %Plan.Par2Rec, at the moment only linearly growing and
                %singletons
                NunidiffofsetRead=1;                
                if length(offsetRead)==1
                    sizeRead=sizRead(1)/2;
                else
                    diffoffsetRead=diff(offsetRead);
                    unidiffoffsetRead=unique(diffoffsetRead)';                                   
                    
                    NunidiffofsetRead=length(unidiffoffsetRead);                    
                    if NunidiffofsetRead~=1
                        %THIS IS A POSSIBILITY IN CASE SOME WEIRD DATA IS INTERLEAVED WITH THE ACTUAL DATA
                        %if NunidiffofsetRead==2
                        %    countdiffoffsetRead=sum(bsxfun(@eq,diffoffsetRead,unidiffoffsetRead),1);
                        %    [~,indUniRelevant]=max(countdiffoffsetRead);
                        %    sizBlock=find(diffoffsetRead==unidiffoffsetRead(setdiff(1:2,indUniRelevant)),1,'first');
                        %    sizeRead=sizBlock*unidiffoffsetRead(indUniRelevant)/2;
                        %    for b=1:sizBlock:length(offsetRead)
                        %        su=fseek(fidraw,offsetRead(b),-1);if su<0;error('Error %s for file %s',ferror(fidraw),rawFile);end
                        %        [aux,su]=fread(fidraw,sizeRead,'int16=>single');if su~=sizeRead;error('Read %d elements for intended %d with error %s for file %s',su,sizeRead,ferror(fidraw),rawFile);end
                        %        rec.(datTyp)=vertcat(rec.(datTyp),aux);
                        %    end                                
                        %else        
                            fclose(fidraw);fprintf('OFFSETS OF TYP %d ARE NOT GROWING LINEARLY FOR FILE %s: SWITCHING TO RECONFRAME READ\n',m,name);rec.(datTyp)=[];rec.Fail=2;recOu{1}=rec;return;
                        %end
                    end           
                    sizeRead=length(offsetRead)*(offsetRead(2)-offsetRead(1))/2;
                end
                if NunidiffofsetRead==1
                    su=fseek(fidraw,offsetRead(1),-1);if su<0;error('Error %s for file %s',ferror(fidraw),rawFile);end
                    [aux,su]=fread(fidraw,sizeRead,'int16=>single');if su~=sizeRead;error('Read %d elements for intended %d with error %s for file %s',su,sizeRead,ferror(fidraw),rawFile);end
                    rec.(datTyp)=vertcat(rec.(datTyp),aux);
                end
            end
        else
            rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==m)=[];
        end
    end
    su=fclose(fidraw);if su<0;error('Error %s for file %s',ferror,rawFile);end
else%USING RECONFRAME
    for m=typ2Rec';datTyp=rec.Plan.Types{m};
        indRead=(rec.Plan.Typ==m);
        if m~=5%We always read all the noise samples if indicated by typ2Rec
            for n=2:rec.Plan.NDims            
                if n~=4 || NC==1;indRead=(indRead & ismember(rec.Plan.Assign{n},rec.Dyn.Ind2Rec{n}));end
            end
        end
        indRead=find(indRead);
        if ~isempty(indRead)
            indReadFull=bsxfun(@plus,NC*(indRead-1)',chancompi);
            rec.Corr.(datTyp)=cell(1,rec.Plan.NCorrs);
            for n=1:rec.Plan.NCorrs
                if n~=2;rec.Corr.(datTyp){n}=rec.Plan.Corr{n}(indRead);else rec.Corr.(datTyp){n}=rec.Plan.Corr{n}(indReadFull(:));end
            end
            rec.Assign.(datTyp)=cell(1,rec.Plan.NDims);
            for n=2:rec.Plan.NDims
                if n~=4 || NC==1;rec.Assign.(datTyp){n}=rec.Plan.Assign{n}(indRead);end
            end                                                
            MR=MRecon(rawFile);
            MR.Parameter.Parameter2Read.typ=uint8(m);
            if m~=5%We always read all the noise samples if indicated by typ2Rec
                for n=2:rec.Plan.NDims;MR.Parameter.Parameter2Read.(rec.Plan.Dims{n})=double(rec.Dyn.Ind2Rec{n});end
            end
            MR.ReadData;      
            %MR.RandomPhaseCorrection;
            %MR.PDACorrection;
            %MR.DcOffsetCorrection;
            %MR.MeasPhaseCorrection;
            %rec.(datTyp)=MR.Data;            
            if ~isempty(MR.Data)
                rec.(datTyp)=[real(MR.Data(:)) imag(MR.Data(:))]';
                rec.(datTyp)=rec.(datTyp)(:);
            end
        else
            rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==m)=[];
        end
    end
    %recOu{1}=rec;
    %return
end    

%CHECK FOR MISSING DATA
if ~isfield(rec,'z') || isempty(rec.z);fprintf('No image data was read so series %s cannot be reconstructed\n',name);recOu{1}=rec;recOu{1}.Fail=1;return;end

typ2Rec=rec.Dyn.Typ2Rec;
%ASSIGN TO STRUCTURE
for m=typ2Rec';datTyp=rec.Plan.Types{m};
    if rec.Dyn.GPU
        if rec.Dyn.GPU~=2;rec.(datTyp)=gpuArray(rec.(datTyp));end
        for n=1:rec.Plan.NCorrs;rec.Corr.(datTyp){n}=gpuArray(rec.Corr.(datTyp){n});end
        for n=1:rec.Plan.NDims;rec.Assign.(datTyp){n}=gpuArray(rec.Assign.(datTyp){n});end        
    end
end

%DETECT SITUATIONS WITH DIFFERENT CHANNEL TYPES AND MODIFY THE STRUCTURE ACCORDINGLY
NCG=length(rec.Par.Labels.CoilNrsPerStack);
assert(NC==1 || NCG==1,'Compressed (%d) channels with non-single (%d) stack of coils',NC,NCG);
%assert(NCG==1 || rec.Dyn.FirstIt,'More than one packets of data for assignment of non-single (%d) stack of coils are not supported (file %s)',NCG,name);
for n=1:NCG;rec.Par.Labels.CoilNrsPerStack{n}=single(rec.Par.Labels.CoilNrsPerStack{n});end
if rec.Dyn.FirstIt && NC==1 && length(rec.Par.Labels.CoilNrsPerStack)>1
    assert(length(rec.Par.Labels.CoilNrsPerStack)<=2,'More than two coil stacks acquired for file %s',name);
    if NCG==2 && length(rec.Par.Labels.CoilNrsPerStack{2})>length(rec.Par.Labels.CoilNrsPerStack{1});surfStack=2;else surfStack=1;end
    %assert(~(length(rec.Par.Labels.CoilNrsPerStack)==2 && length(rec.Par.Labels.CoilNrsPerStack{2})>length(rec.Par.Labels.CoilNrsPerStack{1})),'Second stack having more elements than first stack for file %s',name);
    if isempty(rec.Par.Labels.CoilNrsPerStack);fprintf('Coil Nrs Per Stack is empty\n');rec.Fail=1;return;end
    indSurf=single(rec.Par.Labels.CoilNrsPerStack{surfStack});
    chan=unique(rec.Assign.z{4});

    if length(chan)>16
        rec.Par.Labels.CoilNrsPerStack{1}=find(ismember(chan,indSurf))-1;
        if NCG==2;rec.Par.Labels.CoilNrsPerStack{2}=find(~ismember(chan,indSurf))-1;end
    else
        %THIS HAS BEEN CHANGED TO RECONSTRUCT ANTHONY'S DATASETS. IT MAY BE THAT IT HAS TO BE APPLIED IN ALL CASES BUT HAVEN'T HAD TIME TO TEST IT
        rec.Par.Labels.CoilNrsPerStack{1}=chan(ismember(chan,indSurf));
        if NCG==2;rec.Par.Labels.CoilNrsPerStack{2}=chan(~ismember(chan,indSurf));end
    end    
end

if NCG==0;fprintf('Coil Nrs Per Stack is empty\n');rec.Fail=1;recOu{1}=rec;return;end

%BUILD DIFFERENT STRUCTURES FOR THE DIFFERENT CHANNEL TYPES
recOu=cell(1,NCG);
for n=1:NCG
    recOu{n}=rec;
    recOu{n}.Par.Labels.CoilNrsPerStack=[];    
    recOu{n}.Par.Labels.CoilNrsPerStack{1}=rec.Par.Labels.CoilNrsPerStack{n};
    if n==1
        indS=rec.Assign.z{8}==0;      
        if rec.Par.Mine.Modal~=2;indS(:)=1;end        
        surfaceCoils=unique(gather(rec.Assign.z{4}(indS)));
        if rec.Dyn.FirstIt;recOu{n}.Par.Labels.CoilNrsPerStack{1}=recOu{n}.Par.Labels.CoilNrsPerStack{1}-min(recOu{n}.Par.Labels.CoilNrsPerStack{1});end
        for m=typ2Rec';datTyp=rec.Plan.Types{m};
            indS=rec.Assign.(datTyp){8}==0;
            if rec.Par.Mine.Modal~=2;indS(:)=1;end            
            for l=1:length(surfaceCoils)
                recOu{1}.Assign.(datTyp){4}(indS & rec.Assign.(datTyp){4}==surfaceCoils(l))=recOu{1}.Par.Labels.CoilNrsPerStack{1}(l); 
            end
        end   
    end    
    if n==2 && isempty(rec.Par.Labels.CoilNrsPerStack{2})
        bodyCoils=unique(gather(rec.Assign.z{4}(rec.Assign.z{8}==1)));
        if rec.Dyn.FirstIt;recOu{2}.Par.Labels.CoilNrsPerStack{1}=(0:length(bodyCoils)-1)'+length(rec.Par.Labels.CoilNrsPerStack{1});end
        for m=typ2Rec';datTyp=rec.Plan.Types{m};        
            for l=1:length(bodyCoils)               
                recOu{1}.Assign.(datTyp){4}(rec.Assign.(datTyp){8}==1 & rec.Assign.(datTyp){4}==bodyCoils(l))=recOu{2}.Par.Labels.CoilNrsPerStack{1}(l);       
            end
            recOu{2}.Assign.(datTyp){4}=recOu{1}.Assign.(datTyp){4};
        end       
    end
end

for n=1:NCG
    for m=typ2Rec';datTyp=rec.Plan.Types{m};        
        NP=length(rec.Corr.(datTyp){2});%Number of profiles
        NM=length(rec.(datTyp));
        assert(mod(NM,2*NP)==0,'Number of measures (%d) is not a multiple of number of profiles (%d) in typ %s for file %s',NM,2*NP,datTyp,name);  
                      
        indInt=ismember(recOu{n}.Assign.(datTyp){4},recOu{n}.Par.Labels.CoilNrsPerStack{1});
        if any(indInt~=0)%Otherwise another set of coils may have been used for this calibration data
            recOu{n}.(datTyp)=reshape(rec.(datTyp),[NM/NP NP]);
            recOu{n}.(datTyp)=recOu{n}.(datTyp)(:,indInt);
            recOu{n}.(datTyp)=recOu{n}.(datTyp)(:);
            for l=1:rec.Plan.NCorrs;recOu{n}.Corr.(datTyp){l}=recOu{n}.Corr.(datTyp){l}(indInt);end
            for l=1:rec.Plan.NDims;
                if ~isempty(rec.Assign.(datTyp){l});recOu{n}.Assign.(datTyp){l}=recOu{n}.Assign.(datTyp){l}(indInt);end
            end
        else         
            recOu{n}.(datTyp)=reshape(rec.(datTyp),[NM/NP NP]);
            for l=1:rec.Plan.NCorrs;recOu{n}.Corr.(datTyp){l}=recOu{n}.Corr.(datTyp){l};end
            for l=1:rec.Plan.NDims;recOu{n}.Assign.(datTyp){l}=recOu{n}.Assign.(datTyp){l};end                                                
        end
        %We discard those cases on which the number of coils does not match
        %the total number        
        if length(unique(gather(recOu{n}.Assign.(datTyp){4})))~=length(recOu{n}.Par.Labels.CoilNrsPerStack{1})
            if rec.Dyn.Debug>=1;fprintf('Channel samples (%d) do not correspond to all the possible channels (%d) in typ %s for file %s\n',length(unique(gather(recOu{n}.Assign.(datTyp){4}))),length(recOu{n}.Par.Labels.CoilNrsPerStack{1}),datTyp,name);end
            recOu{n}.Plan.Typ2Rec(recOu{n}.Plan.Typ2Rec==m)=[];
            recOu{n}.Dyn.Typ2Rec(recOu{n}.Dyn.Typ2Rec==m)=[];
            recOu{n}=rmfield(recOu{n},datTyp);
            recOu{n}.Corr=rmfield(recOu{n}.Corr,datTyp);
            recOu{n}.Assign=rmfield(recOu{n}.Assign,datTyp);
        end
    end
end
