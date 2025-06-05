function rec=reconPipeline(rec)

%RECONPIPELINE   Runs the reconstruction pipeline of a given series
%   REC=RECONPIPELINE(REC) runs the reconstruction pipeline going through
%   the different dynamics in the data
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (.Names) and the .lab information (.Par)
%   ** REC is a reconstruction structure with filled reconstruction 
%   information
%

%RECONTRUCTION CONTROL
fprintf('Series %s\n',rec.Names.rawFile);
rec.Fail=0;rec.Anomaly.Detected=zeros(1,2);rec.Anomaly.PercReads=zeros(2,2);rec.Anomaly.PercSlices=zeros(2,2);rec.Anomaly.PercVolumes=zeros(2,2);rec.Anomaly.Coils=[];
rec.Par.Mine.Proce=0;rec=reconPlanner(rec);rec=reconAlgorithm(rec);rec=reconSpecific(rec);
if isfield(rec,'recV')
    for n=1:length(rec.recV)
        rec.recV{n}.Fail=0;rec.recV{n}.Anomaly.Detected=[0 0];rec.Anomaly.PercReads=zeros(2,2);rec.Anomaly.PercSlices=zeros(2,2);rec.Anomaly.PercVolumes=zeros(2,2);rec.Anomaly.Coils=[];
        rec.recV{n}=reconPlanner(rec.recV{n});rec.recV{n}=reconAlgorithm(rec.recV{n});rec.recV{n}=reconSpecific(rec.recV{n}); 
    end
end

if rec.Alg.OnlyJSON;rec.Dyn.Typ2Wri(:)=0;rec.Dyn.Typ2Rec=[];end
read=1;

%RECONSTRUCTION CHUNKS
if ~rec.Fail && ~rec.Alg.OnlyJSON
    Par2Rec=rec.Plan.Par2Rec;Block2Rec=rec.Dyn.Block2Rec;Size=rec.Dyn.Size;
    for e=1:Block2Rec(7):Size(7);rec.Dyn.Ind2Rec{7}=Par2Rec{7}(e:min(e+Block2Rec(7)-1,Size(7)));rec.Dyn.CurItDim(7)=1;%Echoes out of the loop-not sure it is a good idea, as they should cover the whole range of values for adequate sorting of data after the reconstruction    
        if rec.Dyn.Ind2Rec{7}(1)==Par2Rec{7}(1);rec.Dyn.FirstItDim(7)=1;end;
        if rec.Dyn.Ind2Rec{7}(end)==Par2Rec{7}(Size(7));rec.Dyn.LastItDim(7)=1;end
        for b=1:Block2Rec(9):Size(9);rec.Dyn.Ind2Rec{9}=Par2Rec{9}(b:min(b+Block2Rec(9)-1,Size(9)));rec.Dyn.CurItDim(9)=1;%Mixes out of the loop-not sure it is a good idea, as they should cover the whole range of values for adequate sorting of data after the reconstruction        
            if rec.Dyn.Ind2Rec{9}(1)==Par2Rec{9}(1);rec.Dyn.FirstItDim(9)=1;end;
            if rec.Dyn.Ind2Rec{9}(end)==Par2Rec{9}(Size(9));rec.Dyn.LastItDim(9)=1;end
            for c=1:Block2Rec(6):Size(6);rec.Dyn.Ind2Rec{6}=Par2Rec{6}(c:min(c+Block2Rec(6)-1,Size(6)));rec.Dyn.CurItDim(6)=1;
                if rec.Dyn.Ind2Rec{6}(1)==Par2Rec{6}(1);rec.Dyn.FirstItDim(6)=1;end;
                if rec.Dyn.Ind2Rec{6}(end)==Par2Rec{6}(Size(6));rec.Dyn.LastItDim(6)=1;end
                for m=1:Block2Rec(5):Size(5);rec.Dyn.Ind2Rec{5}=Par2Rec{5}(m:min(m+Block2Rec(5)-1,Size(5)));rec.Dyn.CurItDim(5)=1;
                    if rec.Dyn.Ind2Rec{5}(1)==Par2Rec{5}(1);rec.Dyn.FirstItDim(5)=1;end;
                    if rec.Dyn.Ind2Rec{5}(end)==Par2Rec{5}(Size(5));rec.Dyn.LastItDim(5)=1;end
                    for d=1:Block2Rec(11):Size(11);rec.Dyn.Ind2Rec{11}=Par2Rec{11}(d:min(d+Block2Rec(11)-1,Size(11)));rec.Dyn.CurItDim(11)=1;
                        if rec.Dyn.Ind2Rec{11}(1)==Par2Rec{11}(1);rec.Dyn.FirstItDim(11)=1;end;
                        if rec.Dyn.Ind2Rec{11}(end)==Par2Rec{11}(Size(11));rec.Dyn.LastItDim(11)=1;end
                            if all(rec.Dyn.LastItDim);rec.Dyn.LastIt=1;end
                            rec.Dyn.Batch=rec.Dyn.Batch+1;

                            %recGT=reconGyroTools(rec,'K2IP');%To test recons      
                            %rec=recGT;                    
                            %return                            
                            if ismember(rec.Alg.DetectAnomalies,[3 4 6 7]) && rec.Fail==2;rec.Fail=0;end%In this case we will explore the whole cohort                            
                            while read
                                if rec.Dyn.Debug>=1 && rec.Fail==0;fprintf('Reading %s\n',rec.Names.rawFile);tsta=tic;end
                                rec.Read=1;
                                if read==2
                                    if rec.Par.Mine.Modal==10
                                        rec.Dyn.Ind2Rec{11}(end)=[];
                                        if isempty(rec.Dyn.Ind2Rec{11});rec.Fail=1;read=0;end
                                    end
                                    if rec.Par.Mine.Modal==9
                                        rec.Dyn.Ind2Rec{5}(end)=[];
                                        if isempty(rec.Dyn.Ind2Rec{5});rec.Fail=1;read=0;end
                                    end                                    
                                end
                                if read                                    
                                    if ~isfield(rec,'recV')
                                        rec=readData(rec);  
                                        if rec{1}.Fail==2%SWITCHING TO RECONFRAME READ, MAINLY FOR CARDIAC
                                            rawFile=rec{1}.Names.rawFile;
                                            MR=MRecon(rawFile);
                                            MR.ReadData;                                             
                                            MR.RandomPhaseCorrection;
                                            MR.RemoveOversampling;
                                            MR.PDACorrection;
                                            MR.DcOffsetCorrection;
                                            MR.MeasPhaseCorrection;
                                            MR.SortData;
                                            MR.GridData;
                                            rec{1}.Fail=0;
                                            rec=invertData(rec,1);
                                            rec.y=single(MR.Data{1});
                                            rec.N=single(MR.Data{5});
                                            if rec.Dyn.GPU==1;rec.y=gpuArray(rec.y);rec.N=gpuArray(rec.N);end                                           
                                            %NOISE STANDARDIZATION
                                            if (rec.Alg.NoiseStand && rec.Par.Mine.Modal~=2) || rec.Alg.DetectAnomalies;rec.y=standardizeCoils(rec.y,rec.N);end
                                            if any(isnan(rec.y(:)));fprintf('Nan after noise standardization\n');rec.Fail=1;return;end
                                            MR.Data{1}=gather(rec.y);                                            
                                            %MR.ZeroFill;                                            
                                            MR.K2IM;
                                            MR.EPIPhaseCorrection;
                                            MR.K2IP;                                            
                                            MR.GridderNormalization;
                                            rec.y=single(MR.Data{1});
                                            if rec.Dyn.GPU==1;rec.y=gpuArray(rec.y);end
                                            rec.y=ifftshift(rec.y,2);
                                            N=size(rec.y);N(end+1:8)=1;
                                            if N(8)>1 && N(3)==1;perm=1:rec.Plan.NDims;perm([3 8])=[8 3];rec.y=permute(rec.y,perm);end
                                            %MULTIBAND ENCODING
                                            AdHocArray=rec.Par.Mine.AdHocArray;
                                            if isempty(AdHocArray) || all(AdHocArray==0:127);rec.Enc.MB=1;return;end
                                            assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
                                            if AdHocArray(1)==101;rec.Enc.MB=AdHocArray(4);else rec.Enc.MB=1;end
                                            if rec.Enc.MB==1 && AdHocArray(7)<=1 && ((rec.Par.Mine.Modal==9 && ~strcmp(rec.Par.Scan.Technique,'SEEPI')) || rec.Par.Mine.Modal==10)
                                                if isfield(rec.Corr,'P') && size(rec.Corr.P{2},2)==2;rec.Dyn.Typ2Wri(3)=1;end
                                            else
                                                %Z-BLIP ENCODING STRENGTHS
                                                if isempty(rec.Par.Encoding.KzRange)
                                                    rec.Enc.RecSize(3)=size(rec.y,3);
                                                    rec.Enc.FOVSize(3)=size(rec.y,3)*rec.Enc.MB;
                                                    rec.Enc.SlicDist=rec.Par.Labels.VoxelSizes(3)+rec.Par.Labels.SliceGaps(1);%Slice separation
                                                    if AdHocArray(3)==3%TO DO---USE A GENERATEGRID FOR THIS
                                                        sliceGap=AdHocArray(2)/rec.Enc.SlicDist;
                                                        rec.Enc.ExcGrid=(0:rec.Enc.RecSize(3)-1)';
                                                        ExcGridOffs=(0:rec.Enc.MB-1)*sliceGap;
                                                        ExcGridOffs=ExcGridOffs-(rec.Enc.MB-1)*(sliceGap-rec.Enc.RecSize(3))/2;
                                                        rec.Enc.ExcGrid=bsxfun(@plus,rec.Enc.ExcGrid,ExcGridOffs);
                                                        rec.Enc.ExcGrid=reshape(rec.Enc.ExcGrid,[1 rec.Enc.FOVSize(3)]);
                                                    else
                                                        rec.Enc.ExcGrid=0:rec.Enc.FOVSize(3)-1;
                                                    end
                                                end
                                            end            
                                            if rec.Dyn.Debug>=1 && rec.Fail==0;fprintf('Reconstructing %s\n',rec.Names.rawFile);tsta=tic;end
                                            rec=reconstructData(rec);
                                            if rec.Dyn.Debug>=1 && rec.Fail==0;tend=toc(tsta);fprintf('Time reconstructing: %.3f s\n',tend);end

                                            if rec.Dyn.Debug>=1 && rec.Fail==0;fprintf('Writing %s\n',rec.Names.rawFile);tsta=tic;end
                                            rec=writeData(rec);
                                            if rec.Dyn.Debug>=1 && rec.Fail==0;tend=toc(tsta);fprintf('Time writing: %.3f s\n\n',tend);end
                                            
                                            return
                                        end                                        
                                        if isfield(rec{1}.Par.Mine,'InitVol')
                                            rec{1}.Assign.z{11}=rec{1}.Assign.z{11}+rec{1}.Par.Mine.InitVol;
                                            rec{1}.Assign.P{11}=rec{1}.Assign.P{11}+rec{1}.Par.Mine.InitVol;
                                            rec{1}.Dyn.Ind2Rec{11}=rec{1}.Dyn.Ind2Rec{11}+rec{1}.Par.Mine.InitVol;                                            
                                        end                                        
                                    else                            
                                        vVol=rec.Par.Mine.splitInfo(rec.Dyn.Ind2Rec{11}+1);
                                        vVolU=unique(vVol);
                                        for r=vVolU'
                                            rec.recV{r}.Dyn.Ind2Rec=rec.Dyn.Ind2Rec;
                                            rec.recV{r}.Dyn.LastItDim=rec.Dyn.LastItDim;
                                            rec.recV{r}.Dyn.FirstItDim=rec.Dyn.FirstItDim;
                                            rec.recV{r}.Dyn.CurItDim=rec.Dyn.CurItDim;
                                            rec.recV{r}.Dyn.Ind2Rec{11}=rec.Dyn.Ind2Rec{11}(vVol==r)-rec.recV{r}.Par.Mine.InitVol;
                                            rec.recV{r}.Fail=rec.Fail;
                                            recOu=readData(rec.recV{r});
                                            rec=assignReadInfo(recOu{1},rec);
                                            if ~rec.Fail
                                                recOu{1}.Assign.z{11}=recOu{1}.Assign.z{11}+rec.recV{r}.Par.Mine.InitVol;
                                            
                                                typ2Rec=rec.Dyn.Typ2Rec;                                   
                                                for s=typ2Rec';datTyp=rec.Plan.Types{s};
                                                    if r==vVolU(1)
                                                        rec.(datTyp)=recOu{1}.(datTyp);                                            
                                                        for t=1:rec.Plan.NCorrs;rec.Corr.(datTyp){t}=recOu{1}.Corr.(datTyp){t};end
                                                        for t=1:rec.Plan.NDims;rec.Assign.(datTyp){t}=recOu{1}.Assign.(datTyp){t};end
                                                    else                                            
                                                        if ~isempty(recOu{1}.(datTyp));rec.(datTyp)=vertcat(rec.(datTyp),recOu{1}.(datTyp));end
                                                        for t=1:rec.Plan.NCorrs
                                                            if isempty(rec.Corr.(datTyp){t});rec.Corr.(datTyp){t}=recOu{1}.Corr.(datTyp){t};
                                                            elseif ~isempty(recOu{1}.Corr.(datTyp){t});rec.Corr.(datTyp){t}=vertcat(rec.Corr.(datTyp){t},recOu{1}.Corr.(datTyp){t});
                                                            end
                                                        end
                                                        for t=1:rec.Plan.NDims
                                                            if isempty(rec.Assign.(datTyp){t});rec.Assign.(datTyp){t}=recOu{1}.Assign.(datTyp){t};
                                                            elseif ~isempty(recOu{1}.Assign.(datTyp){t});rec.Assign.(datTyp){t}=vertcat(rec.Assign.(datTyp){t},recOu{1}.Assign.(datTyp){t});
                                                            end                                                
                                                        end
                                                    end                                            
                                                end
                                            end
                                        end
                                        for r=1:length(rec.recV)
                                            rec.recV{r}.Dyn.FirstIt=0;
                                            rec.recV{r}=assignReadInfo(rec,rec.recV{r});
                                        end
                                        recOu{1}=rec;
                                        rec=recOu;recOu=[];
                                    end
                                    if rec{1}.Dyn.Debug>=1 && rec{1}.Fail==0;tend=toc(tsta);fprintf('Time reading: %.3f s\n',tend);end

                                    %if isempty(strfind(rec{1}.Names.pathIn,'raw-ingenia'))
                                        if rec{1}.Dyn.Debug>=1 && rec{1}.Fail==0;fprintf('Sorting %s\n',rec{1}.Names.rawFile);tsta=tic;end
                                        rec=sortData(rec);
                                        if rec{1}.Dyn.Debug>=1 && rec{1}.Fail==0;tend=toc(tsta);fprintf('Time sorting: %.3f s\n',tend);end
                                    %end
                                    
                                    if rec{1}.Read==1;read=0;end
                                    if rec{1}.Read==2;rec=rec{1};read=2;end
                                else
                                    recOu=rec;rec=[];
                                    rec{1}=recOu;recOu=[];                                    
                                end
                            end                            
                            read=1;rec{1}.Read=read;

                            if rec{1}.Dyn.Debug>=1 && rec{1}.Fail==0;fprintf('Inverting %s\n',rec{1}.Names.rawFile);tsta=tic;end
                            rec=invertData(rec);
                            if rec.Dyn.Debug>=1 && rec.Fail==0;tend=toc(tsta);fprintf('Time inverting: %.3f s\n',tend);end                            
                            if ismember(rec.Alg.DetectAnomalies,[1 2 3 5 6]) && any(rec.Anomaly.Detected);break;end
                            
                            if rec.Dyn.Debug>=1 && rec.Fail==0;fprintf('Reconstructing %s\n',rec.Names.rawFile);tsta=tic;end
                            rec=reconstructData(rec);
                            if rec.Dyn.Debug>=1 && rec.Fail==0;tend=toc(tsta);fprintf('Time reconstructing: %.3f s\n',tend);end                            
                            rec.Dyn.FirstIt=0;rec.Dyn.CurItDim(:)=0;
                    end;rec.Dyn.FirstItDim(11)=0;rec.Dyn.LastItDim(11)=0;
                end;rec.Dyn.FirstItDim(5)=0;rec.Dyn.LastItDim(5)=0;
            end;rec.Dyn.FirstItDim(6)=0;rec.Dyn.LastItDim(6)=0;
        end;rec.Dyn.FirstItDim(9)=0;rec.Dyn.LastItDim(9)=0;
    end;rec.Dyn.FirstItDim(7)=0;rec.Dyn.LastItDim(7)=0;
end

if ~rec.Fail || (isfield(rec,'x') && ~isempty(rec.x))    
    %ESTIMATE B1 FIELD
    if rec.Par.Mine.Modal==4 && isfield(rec,'x') && ~isempty(rec.x);rec=solveB1(rec);end
    
    if rec.Fail && (~isfield(rec,'Recursive') || ~rec.Recursive)
        fprintf('RECONSTRUCTION ERRORS\n');
        rec.Fail=0;
        %rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==1)=[];
        rec.Dyn.Typ2Rec=[];%Previous line was weird
        rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,12);%WE FORCE WRITING OF RECONSTRUCTION DATA     
        %if isfield(rec,'G');rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,9);end%As it may not be complete  
        if isfield(rec,'G');rec=rmfield(rec,'G');end%As it may not be complete
    end
    if ~isfield(rec,'x') || ~isfield(rec,'Enc')
        if ~rec.Alg.DetectAnomalies;fprintf('RECONSTRUCTION FAILED. LAST MESSAGE SHOULD PROVIDE THE REASON\n');end
    else
        if rec.Dyn.Debug>=1;fprintf('Writing %s\n',rec.Names.rawFile);tsta=tic;end
        rec=writeData(rec);
        if rec.Dyn.Debug>=1;tend=toc(tsta);fprintf('Time writing: %.3f s\n\n',tend);end
    end
else
    if ~rec.Alg.DetectAnomalies && (~isfield(rec,'adding') || rec.adding==0);fprintf('RECONSTRUCTION FAILED. LAST MESSAGE SHOULD PROVIDE THE REASON\n');end
end

function y=assignReadInfo(x,y)
    y.Plan.Typ2Rec=x.Plan.Typ2Rec;
    y.Dyn.Typ2Rec=x.Dyn.Typ2Rec;
    y.Fail=x.Fail;
    y.Par.Labels.CoilNrsPerStack=x.Par.Labels.CoilNrsPerStack;
end

end
