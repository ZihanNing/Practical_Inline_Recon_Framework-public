function rec=sortData(rec)

%SORTDATA   Performs initial corrections to the raw data and places it in 
%the grid used for reconstruction
%   REC=SORTDATA(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction 
%   (rec.Corr.(rec.Plan.Types)) and the informaton for sorting 
%   (rec.Assign(rec.Plan.Types))
%   ** REC is a reconstruction structure with sorted data information 
%   (rec.(rec.Plan.Types)) and removed information for correction 
%   (rec.Corr.(rec.Plan.Types)) and sorting (rec.Assign(rec.Plan.Types))
%

NCG=length(rec);
for ncg=1:NCG
    if rec{ncg}.Fail;return;end

    NC=length(rec{ncg}.Corr.z{2})/length(rec{ncg}.Corr.z{1});
    typ2Rec=rec{ncg}.Dyn.Typ2Rec;name=rec{ncg}.Names.Name; 

    for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};
        %CHECKS ON THE DIFFERENT CHANNELS
        Nm=length(rec{ncg}.(datTyp));%Number of measures            
        assert(mod(Nm,2)==0,'The number of measures is not even, so not clear how to assign quadrature rec{ncg}eptors in typ %s for file %s',datTyp,name);
        Nm=Nm/2;%Number of measures
        NP=length(rec{ncg}.Corr.(datTyp){2});%Number of profiles
        chan=unique(gather(rec{ncg}.Assign.(datTyp){4}));
        if NC==1;NChan=length(chan);else NChan=NC;end   

        assert(mod(Nm,NChan)==0,'Number of measures (%d) divided by number of coils (%d) is not an integer in typ %s for file %s',Nm,NChan,datTyp,name);
        assert(mod(NP,NChan)==0,'Number of profiles (%d) divided by number of coils (%d) is not an integer in typ %s for file %s',NP,NChan,datTyp,name);    

        if rec{ncg}.Alg.CheckGain;rec{ncg}.Max(1)=max(abs(rec{ncg}.(datTyp)(:)));end  
%        max(abs(rec{ncg}.(datTyp)(:)))
                   
%         %CHECKING COMPANDING, DOES NOT PRODUCE A DIFFERENCE         
%          if n==1  
%             nobits=12;
%             x=rec{ncg}.(datTyp);
%             rec{ncg}.(datTyp)=rec{ncg}.(datTyp)*11.1111;
%             rec{ncg}.(datTyp)=compand(rec{ncg}.(datTyp),87.6,2^15,'a/compressor');
%             rec{ncg}.(datTyp)=2^(16-nobits)*round(rec{ncg}.(datTyp)/2^(16-nobits));
%             rec{ncg}.(datTyp)=compand(rec{ncg}.(datTyp),87.6,2^15,'a/expander');
%             rec{ncg}.(datTyp)=rec{ncg}.(datTyp)/11.1111;
%             max(abs(rec{ncg}.(datTyp)(:)-x(:)))
%          end
        
        %QUADRATURE RECEPTION TO COMPLEX DATA
        %if isempty(strfind(rec{ncg}.Names.pathIn,'raw-ingenia'))
            rec{ncg}.(datTyp)=reshape(rec{ncg}.(datTyp),[2 Nm]);
            rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),1,1)+1i*dynInd(rec{ncg}.(datTyp),2,1);%2xSingle to Complex  
        %end
        
        if rec{ncg}.Dyn.GPU==2;rec{ncg}.(datTyp)=gpuArray(rec{ncg}.(datTyp));end

        %COMPLEX MEASURES FOR EACH CHANNEL
        Nm=Nm/NChan;
        NP=NP/NChan;                
        
        if NC==1
            for m=1:rec{ncg}.Plan.NCorrs                                
                rec{ncg}.Corr.(datTyp){m}=resPop(rec{ncg}.Corr.(datTyp){m},[1 2],[NChan NP],[2 1]);        
            end
            for m=2:rec{ncg}.Plan.NDims
                rec{ncg}.Assign.(datTyp){m}=resPop(rec{ncg}.Assign.(datTyp){m},[1 2],[NChan NP],[2 1]);
            end
            chanOrd=unique(gather(rec{ncg}.Assign.(datTyp){4}),'rows');
            assert(size(chanOrd,1)==1,'Profiles in %d different channel orders in typ %s for file %s',size(chanOrd,1),datTyp,name);
        else
            rec{ncg}.Corr.(datTyp){2}=resPop(rec{ncg}.Corr.(datTyp){2},[1 2],[NChan NP],[2 1]);
        end
        %COMPLEX MEASURES FOR EACH CHANNEL TO PROFILES FOR EACH CHANNEL
        if mod(Nm,NP)~=0%Recent modification
            fprintf('Number of measures (%d) divided by number of profiles (%d) is not an integer in typ %s for file %s\n',Nm,NP,datTyp,name);rec{ncg}.Fail=1;return;
        end            
        
        %if isempty(strfind(rec{ncg}.Names.pathIn,'raw-ingenia'))                
            rec{ncg}.(datTyp)=reshape(resPop(rec{ncg}.(datTyp),[1 2],[Nm/NP NChan NP],[1 3 2]),[Nm NChan]);    
        %end
        rec{ncg}.Assign.(datTyp){4}=permute(rec{ncg}.Assign.(datTyp){4},[1 3 4 2]);

        Nm=Nm/NP;
        rec{ncg}.(datTyp)=reshape(rec{ncg}.(datTyp),[Nm NP 1 NChan]);
        for m=1:rec{ncg}.Plan.NCorrs
            rec{ncg}.Corr.(datTyp){m}=reshape(rec{ncg}.Corr.(datTyp){m},[1 NP 1 NChan]);
        end        
        for m=2:rec{ncg}.Plan.NDims     
            rec{ncg}.Assign.(datTyp){m}=reshape(rec{ncg}.Assign.(datTyp){m},[1 NP 1 NChan]);
            rec{ncg}.Assign.(datTyp){m}=dynInd(rec{ncg}.Assign.(datTyp){m},1,2+2*(m~=4)); %We only preserve the information of the first sample
        end
        
        %PROBABLY BETTER IN CASE OF ANOMALY DETECTION TO ADD QUANTIZATION NOISE TO MAKE X NOT CATHEGORICAL
        if rec{ncg}.Alg.DetectAnomalies==1
            NA=size(rec{ncg}.(datTyp));NA(end+1:4)=1;NA(4)=1;
            for s=1:size(rec{ncg}.(datTyp),4);rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),s,4,dynInd(rec{ncg}.(datTyp),s,4)+rand(NA,'like',real(rec{ncg}.(datTyp)))+1i*rand(NA,'like',real(rec{ncg}.(datTyp)))-0.5*(1+1i));end%For loop to prevent memory issues...
        end     
        
        %RANDOM PHASE PARAMETERS
        rec{ncg}.Corr.(datTyp){1}=single(exp(-(2*pi)*1i*double(rec{ncg}.Corr.(datTyp){1})/(2^16-1)));       
        %PDA PARAMETERS
        PDACorrMeth=2;%0->No correc{ncg}tion / 1->Using pda_index and PDAFactors / 2-> using pda_fac (default)
        if PDACorrMeth==0;rec{ncg}.Corr.(datTyp){2}(:)=1;elseif PDACorrMeth==1;rec{ncg}.Corr.(datTyp){2}=reshape(rec{ncg}.Par.Labels.PDAFactors(rec{ncg}.Corr.(datTyp){3}+1),size(rec{ncg}.Corr.(datTyp){3}));end%And if Alg.PDACorr==2 we don't do anything to rec{ncg}.Corr.(datTyp){2}
        rec{ncg}.Corr.(datTyp){2}=single(rec{ncg}.Corr.(datTyp){2});
        
        %RANDOM PHASE AND PDA CORRECTIONS TOGETHER
        %if isempty(strfind(rec{ncg}.Names.pathIn,'raw-ingenia'))
            rec{ncg}.(datTyp)=bsxfun(@times,rec{ncg}.(datTyp),bsxfun(@times,rec{ncg}.Corr.(datTyp){1},rec{ncg}.Corr.(datTyp){2}));   
        %end
        for m=[1 3];rec{ncg}.Corr.(datTyp){m}=[];end
        
        %CHECK RECEIVER GAIN
        if rec{ncg}.Alg.CheckGain;rec{ncg}.Max(2)=gather(max(cat(1,abs(real(rec{ncg}.(datTyp)(:))),abs(imag(rec{ncg}.(datTyp)(:))))));rec{ncg}.Fail=1;return;end                   
        
        %ACTUAL ANOMALY DETECTION        
        if rec{ncg}.Alg.DetectAnomalies==1;rec{ncg}.Anomaly.Detected(1)=anomalyDetection(rec{ncg}.z,rec{ncg}.Names.pathOu,rec{ncg}.Names.path,name,rec{ncg}.Alg.DetectAnomalies);rec{ncg}.Fail=1;return;end
    end
    
    %DCOFFSET PARAMETERS
    if any(rec{ncg}.Dyn.Typ2Rec==5) && rec{ncg}.Dyn.FirstIt;rec{ncg}.Corr.DCOffset=mean(rec{ncg}.N,1);end%We assume noise is the last one
    %if isempty(strfind(rec{ncg}.Names.pathIn,'raw-ingenia'))
        for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};    
            %DCOFFSET CORRECTION
            if any(rec{ncg}.Dyn.Typ2Rec==5);rec{ncg}.(datTyp)=bsxfun(@minus,rec{ncg}.(datTyp),rec{ncg}.Corr.DCOffset);elseif rec{ncg}.Dyn.Debug>=1;fprintf('No noise samples have been acquired for file %s\n',name);end
            %MEASUREMENT PHASE PARAMETERS
            rec{ncg}.Corr.(datTyp){4}=exp(-(pi/2)*1i*single(rec{ncg}.Corr.(datTyp){4}));
            %MEASUREMENT PHASE CORRECTION
            rec{ncg}.(datTyp)=bsxfun(@times,rec{ncg}.(datTyp),rec{ncg}.Corr.(datTyp){4});
            rec{ncg}.Corr.(datTyp){4}=[];  
            %%CHECK THERE IS MEANINGFUL DATA
            %if any(isnan(rec{ncg}.z(:)));fprintf('Nan data after correction\n');rec{ncg}.Fail=1;return;end
        end
    %end
    
    %SORTING THE DATA
    %ARRANGEMENTS FOR CALIBRATION/NAVIGATION DATA
    for n=3:4;datTyp=rec{ncg}.Plan.Types{n};
        if any(typ2Rec==n) && isfield(rec{ncg}.Assign,datTyp)
            stop_sort=0;
            if numel(rec{ncg}.Assign.(datTyp){4})~=numel(rec{ncg}.Assign.z{4});fprintf('Not matched coils for calibration (%d) and data (%d) for file %s\n',numel(rec{ncg}.Assign.(datTyp){4}),numel(rec{ncg}.Assign.z{4}),name);end        
            if ~any(ismember(rec{ncg}.Assign.(datTyp){4}(:),rec{ncg}.Assign.z{4}(:)))         
                fprintf('Not any calibration coil ids (%s) are matched with data coil ids (%s) for file %s. Calibration is not used\n',sprintf('%d ',rec{ncg}.Assign.(datTyp){4}(:)'),sprintf('%d ',rec{ncg}.Assign.z{4}(:)'),name);
                rec{ncg}.Assign=rmfield(rec{ncg}.Assign,datTyp);rec{ncg}=rmfield(rec{ncg},datTyp);rec{ncg}.Corr=rmfield(rec{ncg}.Corr,datTyp);typ2Rec(typ2Rec==3)=[];rec{ncg}.Dyn.Typ2Rec(rec{ncg}.Dyn.Typ2Rec==n)=[];           
                stop_sort=1;
            end
            %stop_sort
            if ~stop_sort  
                NP=numel(rec{ncg}.Assign.(datTyp){2});NZ=numel(rec{ncg}.Assign.z{2});
                if NP>NZ;fprintf('More calibration (%d) than data (%d) was acquired for file %s\n',NP,NZ,name);rec{ncg}.Fail=1;return;end
                if mod(NZ,NP)~=0;fprintf('Calibration profiles do not seem to come from a well-defined block of data, remainder %.5f\n',mod(NZ,NP));rec{ncg}.Fail=1;return;end
                matchOrd=1;
                if n==3
                    for m=5:rec{ncg}.Plan.NDims                               
                        if any(rec{ncg}.Assign.(datTyp){m}~=rec{ncg}.Assign.z{m}(1:NP));
                            indMatched=find(rec{ncg}.Assign.(datTyp){m}==rec{ncg}.Assign.z{m}(1:NP));                    
                            %indMatched
                            if isempty(indMatched)
                                matchOrd=0;
                                fprintf('Ordering of calibration not matched with ordering of data in dim %s for file %s\n',rec{ncg}.Plan.Dims{m},name);
                                rec{ncg}.Assign=rmfield(rec{ncg}.Assign,datTyp);rec{ncg}=rmfield(rec{ncg},datTyp);rec{ncg}.Corr=rmfield(rec{ncg}.Corr,datTyp);typ2Rec(typ2Rec==n)=[];rec{ncg}.Dyn.Typ2Rec(rec{ncg}.Dyn.Typ2Rec==n)=[];rec{ncg}.Plan.Typ2Rec(rec{ncg}.Plan.Typ2Rec==n)=[];
                                break
                            else
                                for l=2:rec{ncg}.Plan.NDims
                                    if l~=4;rec{ncg}.Assign.(datTyp){l}=rec{ncg}.Assign.(datTyp){l}(indMatched);end
                                end   
                                for c=[2 5];rec{ncg}.Corr.(datTyp){c}=dynInd(rec{ncg}.Corr.(datTyp){c},indMatched,2);end
                                rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),indMatched,2);
                                NP=numel(rec{ncg}.Assign.(datTyp){2});
                                if mod(NZ,NP)~=0;fprintf('Calibration profiles do not seem to come from a well-defined block of data, remainder %.5f\n',mod(NZ,NP));rec{ncg}.Fail=1;return;end
                            end
                        end
                    end                    
                    if matchOrd;rec{ncg}.Assign.(datTyp){2}=rec{ncg}.Assign.z{2}(1:NP);end%We assume that kys for calibration correspond to the first loop of kys for data but without blips                    
                else
                    rec{ncg}.Assign.(datTyp){2}=rec{ncg}.Assign.z{2}(1:NZ/NP:end);
                end 
                dims2ExtractAux=1:rec{ncg}.Plan.NDims;dims2ExtractAux=setdiff(dims2ExtractAux,2);
            end
        end
    end
    
    %TRYING TO SAVE MEMORY
    if rec{ncg}.Dyn.GPU==2
        for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};
            rec{ncg}.(datTyp)=gather(rec{ncg}.(datTyp));
        end
    end   
    
    %SORTING BASED ON DETECTING PERIODIC PATTERNS IN THE DATA
    for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};
        %DIMENSIONS OF ASSIGNMENT
        NP=length(rec{ncg}.Assign.(datTyp){2});
        if rec{ncg}.Dyn.GPU==2;rec{ncg}.(datTyp)=gpuArray(rec{ncg}.(datTyp));end
        rec{ncg}.(datTyp)=permuteInit(rec{ncg}.(datTyp));   
        for c=[2 5];rec{ncg}.Corr.(datTyp){c}=permuteInit(rec{ncg}.Corr.(datTyp){c});end

        confliInterrupt=1;    
        while confliInterrupt     
            rec{ncg}.Corr.(datTyp){1}=2*zeros(ones(1,rec{ncg}.Plan.NDims));
            freqDims=cell(1,rec{ncg}.Plan.NDims);          
            
            singleDims=single(zeros(1,rec{ncg}.Plan.NDims));
            sortedDims=singleDims;sortedDims([1 4])=1;
            sizeDims=singleDims;
            rangeDims=single(zeros(2,rec{ncg}.Plan.NDims));        
            ind2Sort=single(ones(1,NP));

            Nsort=1;
            dims2Sort=[];confli1Dims=[];confli2Dims=[];        
            for m=1:rec{ncg}.Plan.NDims                        
                if rec{ncg}.Dyn.GPU;rec{ncg}.Assign.(datTyp){m}=gather(rec{ncg}.Assign.(datTyp){m});end
                if ~sortedDims(m) && ~ismember(m,[1 4])
                    updateSortInfo;
                end
            end        
            if ~isempty(confli1Dims) && all(ismember(confli1Dims,rec{ncg}.Plan.DimsOutLoop))%DWI conflict outer dimensions
                if any(confli1Dims==10)%We drop dim 10           
                    l=10;
                    indPres=true(1,numel(rec{ncg}.Assign.(datTyp){l}));
                    unD=unique(rec{ncg}.Assign.(datTyp){l}(:));
                    indPres(rec{ncg}.Assign.(datTyp){l}(:)==unD(1))=0;
                    rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),indPres,rec{ncg}.Plan.NDims+1);                
                    for c=[2 5];rec{ncg}.Corr.(datTyp){c}=dynInd(rec{ncg}.Corr.(datTyp){c},indPres,rec{ncg}.Plan.NDims+1);end
                    ind2Sort=dynInd(ind2Sort,indPres,2);
                    for l=1:rec{ncg}.Plan.NDims
                        if ~sortedDims(l);rec{ncg}.Assign.(datTyp){l}=dynInd(rec{ncg}.Assign.(datTyp){l},indPres,2);end
                    end
                    confli1DimsPrev=confli1Dims;
                    confli1Dims=[];
                    for m=confli1DimsPrev;updateSortInfo;end
                    if ~isempty(confli1Dims) && all(ismember(confli1Dims,rec{ncg}.Plan.DimsOutLoop))%Previous method for DWI conflict outer dimensions                                
                        C=length(confli1Dims);
                        if C>1
                            N=single(ones(1,C));
                            for c=1:C;N(c)=length(rec{ncg}.Plan.Par2Rec{confli1Dims(c)});end 
                            [~,indm]=max(N);m=confli1Dims(indm);            
                            if rec{ncg}.Dyn.Debug>=1;fprintf('Resolving conflict in dims%s with privileged dim %d\n',sprintf(' %d',confli1Dims),m);end
                            indPres=true(1,numel(rec{ncg}.Assign.(datTyp){m}));
                            for l=confli1Dims
                                if l~=m
                                    unD=unique(rec{ncg}.Assign.(datTyp){l}(:));
                                    histD=histc(rec{ncg}.Assign.(datTyp){l}(:),unD);
                                    [~,indM]=max(histD);%THIS MAY PROVOKE ISSUES IF THE DYNBLOCK2REC ARE MADE TOO SHORT                                                                      
                                    indPres(rec{ncg}.Assign.(datTyp){l}(:)~=unD(indM))=0;       
                                    rec{ncg}.Assign.(datTyp){l}=0;
                                    singleDims(l)=1;sortedDims(l)=1;
                                end
                            end
                            rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),indPres,rec{ncg}.Plan.NDims+1);                
                            for c=[2 5];rec{ncg}.Corr.(datTyp){c}=dynInd(rec{ncg}.Corr.(datTyp){c},indPres,rec{ncg}.Plan.NDims+1);end
                            ind2Sort=dynInd(ind2Sort,indPres,2);
                            for l=1:rec{ncg}.Plan.NDims
                                if ~sortedDims(l);rec{ncg}.Assign.(datTyp){l}=dynInd(rec{ncg}.Assign.(datTyp){l},indPres,2);end
                            end
                            confli1Dims=[];
                            updateSortInfo;
                        end
                    end
                end
            end
            if ~isempty(confli1Dims) && all(ismember(confli1Dims,[2 7]))%DWI conflict dual echo or partial calibration data
                if all(ismember(confli1Dims,2))%Conflict for calibration data                  
                    if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict in dim 2 that needs recoding\n');end
                else%DWI conflict dual echo                                        
                    if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict in dual echo data (dims 2 and 7) that needs recoding\n');end
                end  
                confliInterrupt=0;
            elseif ~isempty(confli1Dims) && all(ismember(confli1Dims,[2 3]))%Volumetric scans conflict
                if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict in dims 2 and 3 that needs to be resolved\n');end
                confliInterrupt=0;
            elseif ~isempty(confli1Dims) && all(ismember(confli1Dims,12))
                if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict in dim 12 that needs to be resolved\n');end
                confliInterrupt=0;
            elseif ~isempty(confli1Dims) && all(ismember(confli1Dims,[2 3 12]))%Volumetric scans with interleaved averages conflict
                if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict in dims 2, 3 and 12 that needs to be resolved\n');end
                confliInterrupt=0;
            elseif ~isempty(confli1Dims)            
                if rec{ncg}.Dyn.Debug>=1;fprintf('There is a conflict probably related with interrupted scan:%s\n',sprintf(' %d',confli1Dims));end
                if any(confli1Dims==5)
                    indDisc=find(rec{ncg}.Assign.(datTyp){5}~=max(rec{ncg}.Assign.(datTyp){5}));
                    rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),indDisc,rec{ncg}.Plan.NDims+1);
                    for c=[2 5];rec{ncg}.Corr.(datTyp){c}=dynInd(rec{ncg}.Corr.(datTyp){c},indDisc,rec{ncg}.Plan.NDims+1);end
                    ind2Sort=dynInd(ind2Sort,indDisc,2);
                    NP=numel(ind2Sort);
                    for l=1:rec{ncg}.Plan.NDims
                        if ~sortedDims(l);rec{ncg}.Assign.(datTyp){l}=dynInd(rec{ncg}.Assign.(datTyp){l},indDisc,2);end                
                    end
                elseif rec{ncg}.Read==1
                    fprintf('Unresolved conflict when sorting%s',sprintf(' %d',confli1Dims));
                    if ismember(rec{ncg}.Par.Mine.Modal,9:10)                    
                        fprintf('. Trying to reread\n');rec{ncg}.Read=rec{ncg}.Read+1;return;
                    else
                        fprintf('\n');
                        rec{ncg}.Fail=1;return;
                    end
                else
                    fprintf('Unresolved conflict when sorting%s; not reconstructing\n',sprintf(' %d',confli1Dims));rec{ncg}.Fail=1;return;
                end
            elseif ~isempty(confli2Dims)      
                error('Unresolved conflict when sorting:%s',sprintf(' %d',confli2Dims));
            else
                confliInterrupt=0;
            end
        end     
        
        %IF SEVERAL VOLUMES ARE TO BE RECONSTRUCTED AND THE ACQ WAS
        %INTERRUPTED, THE DATA FROM THE LAST VOLUME IS ALREADY DELETED
        %HERE, OTHERWISE IT WILL BE IN INVERTDATA. IN THE FIRST CASE ERRORS
        %IN THE RECONSTRUCTION WON'T BE LOGGED BUT THEY WILL IN THE SECOND
        ind2Sort=double(ind2Sort);
        [per,rep]=seqperiod(ind2Sort);
        if mod(rep,1)~=0;fprintf('An interrupted sequence of data may have been encountered for dimension %d in typ %s for file %s. The data will not be reconstructed.\n',m,datTyp,name);rec{ncg}.Fail=1;return;end
        rep=sum(diff(ind2Sort)~=0)+1;     
        NP=numel(ind2Sort);            
                               
        if per<NP || rep==Nsort%Inside/Outside the loop: we can compress           
            [~,indS]=sort(ind2Sort);
            outLoop=0;
            if rep==Nsort
                outLoop=1;
                indS=reshape(indS,[NP/Nsort Nsort])';indS=indS(:)';
            end            
            rec{ncg}.(datTyp)=sortAndReshapeAndCompressAndPopulate(rec{ncg}.(datTyp),indS,Nsort,0,sizeDims(dims2Sort),dims2Sort);
            for c=[2 5];rec{ncg}.Corr.(datTyp){c}=sortAndReshapeAndCompressAndPopulate(rec{ncg}.Corr.(datTyp){c},indS,Nsort,0,sizeDims(dims2Sort),dims2Sort);end
            dims2Extract=[];
            for m=1:rec{ncg}.Plan.NDims
                if ~sortedDims(m)
                    rec{ncg}.Assign.(datTyp){m}=permuteInit(rec{ncg}.Assign.(datTyp){m});
                    isdis2Sort=single(ismember(m,dims2Sort));
                                       
                    if ~isdis2Sort                   
                        rec{ncg}.Assign.(datTyp){m}=sortAndReshapeAndCompressAndPopulate(rec{ncg}.Assign.(datTyp){m},indS,Nsort,1,sizeDims(dims2Sort),dims2Sort);%%RECENT CHANGE WITH REGARD TO PREVIOUS ONE!!                                        
                    else
                        rec{ncg}.Assign.(datTyp){m}=sortAndReshapeAndCompressAndPopulate(rec{ncg}.Assign.(datTyp){m},1:length(indS),Nsort,0,sizeDims(dims2Sort),dims2Sort);%%RECENT CHANGE WITH REGARD TO PREVIOUS ONE!!                    
                    end
                    if isdis2Sort
                        dims2ExtractAux=1:rec{ncg}.Plan.NDims;dims2ExtractAux=setdiff(dims2ExtractAux,m);
                        aux=dynInd(rec{ncg}.Assign.(datTyp){m},ones(1,rec{ncg}.Plan.NDims-1),dims2ExtractAux);
                        if all(diff(aux)>=0);dims2Extract=horzcat(dims2Extract,m);end
                    end
                    if isdis2Sort;rec{ncg}.Corr.(datTyp){1}(m)=m;else rec{ncg}.Corr.(datTyp){1}(m)=rec{ncg}.Plan.NDims+2;end
                end
            end
            for m=1:rec{ncg}.Plan.NDims
                if ~sortedDims(m)
                    if ismember(m,dims2Sort)
                        if ismember(m,dims2Extract)
                            rec{ncg}.Assign.(datTyp){m}=0;
                        else
                            rec{ncg}.Assign.(datTyp){m}=dynInd(rec{ncg}.Assign.(datTyp){m},ones(1,length(dims2Extract)),dims2Extract);
                            rec{ncg}.Assign.(datTyp){m}=rec{ncg}.Assign.(datTyp){m}(:)';
                        end
                    end
                end
            end
        else
            fprintf('The extracted dimensions for sorting are neither fully inside or outside the loop\n');rec{ncg}.Fail=1;return;
        end
        %REORDER SAMPLES WITH POTENTIALLY WRONG LABELS (PARTICULARLY IN SAFE, BUT PERHAPS NOT ONLY)
        if n~=5 && ~isempty(confli1Dims) && all(ismember(confli1Dims,2)) %&& rec{ncg}.Par.Mine.SAFE 
            rec{ncg}.Assign.(datTyp){2}=permuteAlte(rec{ncg}.Assign.(datTyp){2});
            rec{ncg}.(datTyp)=permuteAlte(rec{ncg}.(datTyp));            
            for c=[2 5];rec{ncg}.Corr.(datTyp){c}=permuteAlte(rec{ncg}.Corr.(datTyp){c});end
            rec{ncg}.Corr.(datTyp){1}(2)=2;
        end
        if rec{ncg}.Dyn.Debug>=2 && n~=5;printResults;end    
        if rec{ncg}.Dyn.GPU==2;rec{ncg}.(datTyp)=gather(rec{ncg}.(datTyp));end
    end
end

function x=permuteInit(x)
    perm=1:rec{ncg}.Plan.NDims+2;perm(rec{ncg}.Plan.NDims+1)=2;perm(2)=rec{ncg}.Plan.NDims+1;
    x=permute(x,perm);
end

function x=permuteAlte(x)
    perm=1:rec{ncg}.Plan.NDims+2;perm(rec{ncg}.Plan.NDims+2)=2;perm(2)=rec{ncg}.Plan.NDims+2;
    x=permute(x,perm);
end

function updateSortInfo
    if length(rec{ncg}.Assign.(datTyp){m})>1;assUn=unique(rec{ncg}.Assign.(datTyp){m});else assUn=rec{ncg}.Assign.(datTyp){m}(1);end            
    sizeDims(m)=length(assUn);
    if sizeDims(m)==1
        rec{ncg}.Assign.(datTyp){m}=0;
        rec{ncg}.Corr.(datTyp){1}(m)=0;
        singleDims(m)=1;sortedDims(m)=1;
    else
        rangeDims(:,m)=[min(assUn);max(assUn)];
        freqDims{m}=histc(rec{ncg}.Assign.(datTyp){m}(:),rangeDims(1,m):rangeDims(2,m));
        [per,rep]=seqperiod(double(rec{ncg}.Assign.(datTyp){m}(:)));
        conf=0;
        if length(unique(freqDims{m}))~=1 || (~(per==1 || rep==1) && m==12 && rec{ncg}.Par.Mine.Modal==7)%Two averages on 3D data
            confli1Dims=[confli1Dims m];conf=1;
        end
        if sizeDims(m)~=(diff(rangeDims(:,m))+1)
            confli2Dims(m)=[confli2Dims m];conf=1;
        end

        if ~conf                     
            dims2Sort=[dims2Sort m];
            ind2Sort=ind2Sort+(rec{ncg}.Assign.(datTyp){m}-rangeDims(1,m))*Nsort;
            Nsort=Nsort*sizeDims(m);
            if rec{ncg}.Dyn.Debug>=1;fprintf('Well-behaved dim %d in typ %s for file %s\n',m,datTyp,name);end
        end
    end
end

function x=sortAndReshapeAndCompressAndPopulate(x,indS,Nsort,comp,Nd2sort,d2sort)
    N=size(x);N(end+1:4)=1; 
    x=dynInd(x,indS,rec{ncg}.Plan.NDims+1);
    x=resPop(x,rec{ncg}.Plan.NDims+1,[Nsort NP/Nsort],rec{ncg}.Plan.NDims+(1:2));
    x=dynInd(x,1,rec{ncg}.Plan.NDims+comp);
    if comp==0 && ~isempty(d2sort);x=resPop(x,rec{ncg}.Plan.NDims+1,Nd2sort,d2sort);end
end

function printResults
    fprintf('Size of data typ %s:%s\n',datTyp,sprintf(' %d',size(rec{ncg}.(datTyp))));    
    %fprintf('Size of polarity typ %s:%s\n',datTyp,sprintf(' %d',size(rec{ncg}.Corr.(datTyp){5})));            
    fprintf('Dims of assign typ %s:%s\n',datTyp,sprintf(' %d',rec{ncg}.Corr.(datTyp){1}))    
end

end
