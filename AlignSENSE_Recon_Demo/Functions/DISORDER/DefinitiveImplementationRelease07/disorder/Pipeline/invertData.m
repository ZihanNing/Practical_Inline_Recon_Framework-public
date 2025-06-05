function recOu=invertData(rec,metaData)

%INVERTDATA   Performs the basic Fourier inversion of acquired data. This
%may include Gridding and Nyquist ghosting correction. In addition, it
%builds the information required for pseudo-inverse reconstruction
%   RECOU=INVERTDATA(REC,{METADATA})
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction 
%   (rec.Corr.(rec.Plan.Types)) and the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types))
%   *{METADATA} is a flag to only compute the inversion information 
%   without performing the actual inversion, defaults to 0
%   ** RECOU is a reconstruction structure with backtransformed data 
%   (rec.(rec.PlanTypes)) and encoding (rec.Enc) information
%

if nargin<2 || isempty(metaData);metaData=0;end

NCG=length(rec);
for ncg=1:NCG
    Enc=[];
    if ncg==1 && rec{1}.Fail;recOu=rec{1};return;end

    typ2Rec=rec{ncg}.Dyn.Typ2Rec;name=rec{ncg}.Names.Name;gpu=rec{ncg}.Dyn.GPU;
    if ~strcmp(rec{ncg}.Par.Scan.AcqMode,'Cartesian');rec{ncg}.Fail=1;fprintf('File %s is a %s scan, which cannot be reconstructed\n',name,rec{ncg}.Par.Scan.AcqMode);return;end

    %MODIFICATION TO HANDLE MULTI-ECHO DATA (2018_10_25/DE_58031)
    NKYE=size(rec{ncg}.Par.Encoding.KyRange,1);
    if NKYE>2
        rec{ncg}.Par.Encoding.KyRange(3:NKYE,:)=repmat(rec{ncg}.Par.Encoding.KyRange(2,:),[NKYE-2 1]);
        if size(rec{ncg}.Par.Labels.PartialFourierFactors(:,2),1)==NKYE;rec{ncg}.Par.Labels.PartialFourierFactors(3:NKYE,2)=repmat(rec{ncg}.Par.Labels.PartialFourierFactors(2,2),[NKYE-2 1]);end%Modification to handle multi-echo data
    end    
    
    %BUILD THE FOURIER TRANSFORM
    Enc.kRange={rec{ncg}.Par.Encoding.KxRange,rec{ncg}.Par.Encoding.KyRange,rec{ncg}.Par.Encoding.KzRange};%Sampled ranges
    
    Enc.rRange={rec{ncg}.Par.Encoding.XRange,rec{ncg}.Par.Encoding.YRange,rec{ncg}.Par.Encoding.ZRange};%Sampled ranges
    AuxOvers={rec{ncg}.Par.Encoding.KxOversampling,rec{ncg}.Par.Encoding.KyOversampling,rec{ncg}.Par.Encoding.KzOversampling};%Sampled ranges
    for m=2:3
        if isempty(Enc.kRange{m});Enc.kRange{m}=Enc.kRange{1};Enc.kRange{m}(:)=0;end
        if isempty(Enc.rRange{m});Enc.rRange{m}=Enc.rRange{1};Enc.rRange{m}(:)=0;end
        if isempty(AuxOvers{m});AuxOvers{m}=AuxOvers{1};AuxOvers{m}(:)=1;end    
    end
    Enc.KOvers=single(zeros(1,3));
    for m=1:3
        assert(size(unique(AuxOvers{m},'rows'),1)==1,'The (%d) oversampling factors for dimension %d are not common among instances',size(AuxOvers{m},1),m);
        Enc.KOvers(m)=AuxOvers{m}(1);
    end        
    if size(rec{ncg}.Par.Labels.PartialFourierFactors,1)==rec{ncg}.Par.Encoding.NrEchoes
        Enc.PartFourier=rec{ncg}.Par.Labels.PartialFourierFactors;
    elseif size(unique(rec{ncg}.Par.Labels.PartialFourierFactors,'rows'),1)~=1
        error('The (%d) half scan factors are not common among instances',size(rec{ncg}.Par.Labels.PartFourier,1));
    else
        Enc.PartFourier=rec{ncg}.Par.Labels.PartialFourierFactors(1,:);
    end
    assert(size(unique(rec{ncg}.Par.Labels.SENSEFactor,'rows'),1)==1,'The SENSE factors change in different stacks or other structures for file %s',name);
    rec{ncg}.Par.Labels.SENSEFactor=rec{ncg}.Par.Labels.SENSEFactor(1,:);
    Enc.SamRate=Enc.KOvers./rec{ncg}.Par.Labels.SENSEFactor;    
        
    for m=1:3%TODO---THERE SEEMS TO BE A RESIDUAL SUBPIXEL SHIFT AMONG MODALITIES, BUT THIS MATCHES SCANNER RECONSTRUCTIONS. WE MAY WANT TO IMPROVE OVER THIS          
        if size(unique(rec{ncg}.Par.Labels.PartialFourierFactors(:,m)),1)==size(Enc.kRange{m},1)
            Enc.AcqSize(m)=round(max(bsxfun(@rdivide,diff(Enc.kRange{m},1,2)+1,unique(rec{ncg}.Par.Labels.PartialFourierFactors(:,m))),[],1));%This tells the number of points to include in the oversampled/undersampled spectral grid
        else%Introduced for some datasets on '2018_03_27/FE_330'
            if numDims(Enc.kRange{m})==2 && size(unique(Enc.kRange{m},'rows'),1)==1%Modification for some datasets on '2018_05_14/MA_11130'
                Enc.AcqSize(m)=round(max(bsxfun(@rdivide,diff(dynInd(Enc.kRange{m},1,1),1,2)+1,rec{ncg}.Par.Labels.PartialFourierFactors(:,m)),[],1));%This tells the number of points to include in the oversampled/undersampled spectral grid
            else
                Enc.AcqSize(m)=round(max(bsxfun(@rdivide,diff(Enc.kRange{m},1,2)+1,rec{ncg}.Par.Labels.PartialFourierFactors(:,m)),[],1));%This tells the number of points to include in the oversampled/undersampled spectral grid
            end                       
        end
    end
    %Enc.KOvers
    %Enc.SamRate   
    %rec{ncg}.Par.Encoding.KzOversampling
    %pause
    Enc.SamRate=Enc.SamRate./abs(rec{ncg}.Alg.OverDec);
    Enc.RecSize=round(Enc.AcqSize./max(Enc.SamRate,1)-1e-3);%In case there is oversampling, we crop the data (-1e-3 to prevent numerical instabilities) 
    Enc.FOVSize=round(Enc.RecSize./min(Enc.SamRate,1)-1e-3);%In case there is subsampling we extend the FOV (-1e-3 to prevent numerical instabilities)
    Enc.NDimsGR=2+(Enc.RecSize(3)~=1);
    Enc.Accel=prod(Enc.FOVSize./Enc.RecSize);
    Enc.AcqVoxelSize=min(rec{ncg}.Par.Scan.AcqVoxelSize,[],1);
    
    if strcmp(rec{ncg}.Par.Labels.FastImagingMode,'EPI')
        Enc.MPSOffcentres=bsxfun(@times,rec{ncg}.Par.Scan.MPSOffcentresMM,1./rec{ncg}.Par.Scan.AcqVoxelSize);
    else
        for m=1:3;Enc.MPSOffcentres(:,m,:)=bsxfun(@times,sum(Enc.rRange{m},2)/2,bsxfun(@times,rec{ncg}.Par.Scan.RecVoxelSize(:,m),1./rec{ncg}.Par.Scan.AcqVoxelSize(:,m)));end
    end
    if size(unique(dynInd(Enc.MPSOffcentres,1,3),'rows'),1)~=1;fprintf('The offcentres change in different stacks or other structures for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
    Enc.MPSOffcentres=dynInd(Enc.MPSOffcentres,1,1);
    
    NPE=size(Enc.MPSOffcentres,3);        
    for p=1:NPE
        Enc.rGrid{p}=generateGrid(Enc.RecSize,gpu,1./max(Enc.SamRate,1),[],-dynInd(Enc.MPSOffcentres,p,3)); 
        Enc.rGridAcq{p}=generateGrid(Enc.AcqSize,gpu,1,[],-dynInd(Enc.MPSOffcentres,p,3));
        %rec{ncg}.Par.Scan.MPSOffcentres(1)
        
        if ~isempty(rec{ncg}.Par.Labels.NusEncNrs);Enc.kGrid{1}=rec{ncg}.Par.Labels.NusEncNrs';else Enc.kGrid{1}=min(Enc.kRange{1}(:,1)):max(Enc.kRange{1}(:,2));end
        for m=2:3;Enc.kGrid{m}=min(Enc.kRange{m}(:,1)):max(Enc.kRange{m}(:,2));end
        dkGrid=cell(1,3);dkGridAcq=cell(1,3);
        for m=1:3
            if gpu;Enc.kGrid{m}=gpuArray(Enc.kGrid{m});end    
            [~,dkGrid{m}]=sincKernel(Enc.kGrid{m}(:)/max(Enc.SamRate(m),1),Enc.kGrid{m}(:)/max(Enc.SamRate(m),1),gpu);                           
            [Enc.DFTM{p}{m},Enc.DFTMH{p}{m}]=buildDFTM(Enc.rGrid{p}{m}(:),Enc.kGrid{m}(:),[],dkGrid{m}(:));
            [~,dkGridAcq{m}]=sincKernel(Enc.kGrid{m}(:),Enc.kGrid{m}(:),gpu);              
            [Enc.DFTMAcq{p}{m},Enc.DFTMHAcq{p}{m}]=buildDFTM(Enc.rGridAcq{p}{m}(:),Enc.kGrid{m}(:),[],dkGridAcq{m}(:));
        end
    end
    %FILL KRANGES FOR DENOSING
    if rec{ncg}.Par.Mine.Proce==1 && ~isempty(rec{ncg}.Par.Mine.Signs)
        NPE=size(rec{ncg}.Par.Mine.Signs,3);   
        for n=1:2;kRange(n,:)=Enc.kRange{n}(1,:);end
        for n=1:3
            assert(size(Enc.kRange{n},1)==1,'%d alternating encoding structures not using external file with %d encoding structures',size(Enc.kRange{n},1),NPE);
            Enc.kRange{n}=repmat(Enc.kRange{n},[NPE 1]);
        end
        kRange=repmat(kRange,[1 1 NPE]);
        for n=1:NPE
            perm=(1:2)';
            perm(1:2)=rec{ncg}.Par.Mine.Signs(:,:,n)*perm(1:2);
            kRange(:,:,n)=kRange(abs(perm),:,n);
            for m=1:2
                if perm(m)<0
                    kRange(m,:,n)=-kRange(m,:,n);
                    if mod(Enc.AcqSize(abs(perm(m))),2)==0;kRange(m,:,n)=kRange(m,:,n)-1;end
                end
            end            
        end        
        kRange=sort(kRange,2);
        for n=1:2;Enc.kRange{n}=permute(kRange(n,:,:),[3 2 1]);end
    end

    if metaData
        rec{ncg}.Enc=Enc;recOu=rec{ncg};return;
    %elseif metaData==2
    %    
    end

    %CHECKS GHOSTING CORRECTION
    nDirs=unique(rec{ncg}.Corr.z{5}(:))';
    assert(length(nDirs)<=2,'More than 2 (%d) polarities in the data for file %s',length(nDirs),name);
    assert(all(ismember(nDirs,[1 -1])),'Polarities are not given as pair [-1,1] but as (%s) for file %s',sprintf('%.2f ',nDirs),name);

    %MULTIPLE PHASE ENCODES
    voInds=rec{ncg}.Dyn.Ind2Rec{11}+1;%Volumetric indexes
    indPEV=cell(1,NPE);
    for p=1:NPE
        %if ~isempty(rec{ncg}.Par.Mine.diInfo);indPEV{p}=find(rec{ncg}.Par.Mine.diInfo(voInds,5)==rec{ncg}.Par.Mine.pedsUn(p));else indPEV{p}=1:length(voInds);end
        if ~isempty(rec{ncg}.Par.Mine.diInfo) && any(voInds>size(rec{ncg}.Par.Mine.diInfo,1));fprintf('More volumes in the data than in the direction file for %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
        if ~isempty(rec{ncg}.Par.Mine.diInfo);indPEV{p}=find(rec{ncg}.Par.Mine.diInfo(voInds,5)==p);else indPEV{p}=1:length(voInds);end
    end   

    %TRANSFORM BACK THE READOUT
    if rec{ncg}.Dyn.GPU==2;rec{ncg}.N=gpuArray(rec{ncg}.N);end%From saving memory
    for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};
        if n~=5
            if rec{ncg}.Dyn.GPU==2;rec{ncg}.(datTyp)=gpuArray(rec{ncg}.(datTyp));end%From saving memory
            dims=1:rec{ncg}.Plan.NDims+2;
            if rec{ncg}.Corr.(datTyp){1}(2)==0 || numel(rec{ncg}.Corr.(datTyp){5})/size(rec{ncg}.Corr.(datTyp){5},4)==1
                if rec{ncg}.Dyn.Debug>=1;fprintf('Unclear interpretation of calibration as only a single PE is present in typ %s for file %s\n',datTyp,name);end
                rec{ncg}.Corr.(datTyp){1}(2)=2;
            end
            
            %SAVE RAW BEFORE ANOMALY DETECTION
            if n==1 && ncg==1 && rec{ncg}.Alg.SaveRaw==5;writeRaw(rec{ncg});end
            
            %ANOMALY DETECTION
            if rec{ncg}.Alg.DetectAnomalies>=2 && n==1
                if ~isfield(rec{ncg},'maxProfRisk');rec{ncg}.maxProfRisk=[];end
                for a=1:size(rec{ncg}.z,11)
                    for b=1:size(rec{ncg}.z,7)
                        for c=1:size(rec{ncg}.z,5)
                            if rec{ncg}.Alg.DetectAnomalies<5
                                [rec{ncg}.Anomaly.Detected(1),rec{ncg}.maxProfRisk(end+1)]=anomalyDetection(dynInd(rec{ncg}.z,[c b a],[5 7 11]),rec{ncg}.Names.pathOu,rec{ncg}.Names.path,name,rec{ncg}.Alg.DetectAnomalies,rec{ncg}.N,dynInd(rec{ncg}.Corr.z{2},[c b a],[5 7 11]),dynInd(rec{ncg}.Corr.z{5},[c b a],[5 7 11]),rec{ncg}.Assign,rec{ncg}.Anomaly.Detected(1),[rec{ncg}.Dyn.Ind2Rec{5}(c) rec{ncg}.Dyn.Ind2Rec{7}(b) rec{ncg}.Dyn.Ind2Rec{11}(a)]);
                            else
                                [rec{ncg}.Anomaly.Detected(2),rec{ncg}.Anomaly.PercReads(:,2),rec{ncg}.Anomaly.PercSlices(:,2),rec{ncg}.Anomaly.PercVolumes(:,2),rec{ncg}.Anomaly.Coils]=...
                                    detuningDetection(dynInd(rec{ncg}.z,[c b a],[5 7 11]),rec{ncg}.Names.pathOu,rec{ncg}.Names.path,name,rec{ncg}.Alg.DetectAnomalies,rec{ncg}.N,...
                                    dynInd(rec{ncg}.Corr.z{2},[c b a],[5 7 11]),dynInd(rec{ncg}.Corr.z{5},[c b a],[5 7 11]),rec{ncg}.Assign,...
                                    rec{ncg}.Anomaly.Detected(2),rec{ncg}.Anomaly.PercReads(:,2),rec{ncg}.Anomaly.PercSlices(:,2),rec{ncg}.Anomaly.PercVolumes(:,2),rec{ncg}.Anomaly.Coils,...
                                    [rec{ncg}.Dyn.Ind2Rec{5}(c) rec{ncg}.Dyn.Ind2Rec{7}(b) rec{ncg}.Dyn.Ind2Rec{11}(a)]);
                            end                                
                            if ismember(rec{ncg}.Alg.DetectAnomalies,[2 3 5 6]) && any(rec{ncg}.Anomaly.Detected>0);break;end
                        end
                        if ismember(rec{ncg}.Alg.DetectAnomalies,[2 3 5 6]) && any(rec{ncg}.Anomaly.Detected>0);break;end
                    end
                    if ismember(rec{ncg}.Alg.DetectAnomalies,[2 3 5 6]) && any(rec{ncg}.Anomaly.Detected>0);break;end
                end
                rec{ncg}.Fail=2;recOu=rec{ncg};return;
            end

            %NOISE STANDARDIZATION
            if (rec{ncg}.Alg.NoiseStand && rec{ncg}.Par.Mine.Modal~=2) || rec{ncg}.Alg.DetectAnomalies;rec{ncg}.(datTyp)=standardizeCoils(rec{ncg}.(datTyp),rec{ncg}.N);end                    
            if any(isnan(rec{ncg}.(datTyp)(:)));fprintf('Nan after noise standardization\n');rec{ncg}.Fail=1;recOu=rec{ncg};return;end

            %TRANSFORM BACK THE READOUT            
            aux=gather(mapMat(rec{ncg}.Corr.(datTyp){5},rec{ncg}.Corr.(datTyp){1}(2)));
            assert(size(unique(aux,'rows'),1)<=2,'Polarity information of calibration not consistent among different EPI samples in typ %s for file %s',datTyp,name);aux=[]; 

            dimsInt=setdiff(dims,rec{ncg}.Corr.(datTyp){1}([2 7]));
            rec{ncg}.Corr.(datTyp){5}=dynInd(rec{ncg}.Corr.(datTyp){5},ones(1,length(dimsInt)),dimsInt);

            if ~(rec{ncg}.Corr.(datTyp){1}(2)<7);E=0;else E=rec{ncg}.Plan.Par2Rec{7}';end            
            if NPE==1 || n~=4
                if n==4;E=0:size(rec{ncg}.(datTyp),7)-1;end
                for e=E%Different echoes
                    ex=ones(1,rec{ncg}.Plan.NDims+1);ex(7-single((rec{ncg}.Corr.(datTyp){1}(2))<7))=e+1;   
                    aux=dynInd(rec{ncg}.Corr.(datTyp){5},ex,setdiff(dims,rec{ncg}.Corr.(datTyp){1}(2)));                             
                    for nd=nDirs
                        for p=1:NPE
                            if ismember(n,[1 4]);indV=indPEV{p};elseif p==rec{ncg}.Par.Mine.Nat;indV=1;else indV=[];end
                            if n==4;indV=indV(1:min(length(indV),size(rec{ncg}.(datTyp),11)));end
                            if ~isempty(indV)
                                x=dynInd(rec{ncg}.(datTyp),{aux(:)==nd,e+1,indV},[rec{ncg}.Corr.(datTyp){1}(2) 7 11]);
                                if ~isempty(x)
                                    N=size(x);
                                    N(1)=Enc.RecSize(1);

                                    if rec{ncg}.Alg.PlugNoise==1 && n==1;x=plugNoise(x);end
                                    if nd>0;A=Enc.DFTMH{p}{1};else A=conj(Enc.DFTMH{p}{1});end     
                                    x=aplGPU(A,x,1)/sqrt(trace(A*A')/size(A,1));
                                    if rec{ncg}.Alg.PlugNoise==2 && n==1;x=plugNoise(x);end
                                    rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),{1:N(1),aux(:)==nd,e+1,indV},[1 rec{ncg}.Corr.(datTyp){1}(2) 7 11],x);
                                end
                            end
                        end
                    end
                end
            end
            x=[];rec{ncg}.(datTyp)=dynInd(rec{ncg}.(datTyp),1:N(1),1);
        end
    end
    
    %%NAVIGATOR-BASED CORRECTIONS OF MULTI-SHOT DWI
    %%%THIS WORKS WELL FOR PHANTOM DATA BUT NOT FOR IN VIVO DATA. DON'T UNDERSTAND WHY
    if isfield(rec{ncg},'C') && rec{ncg}.Alg.PhasCorRef && size(rec{ncg}.z,7)==1%Only for single echo data
        rec{ncg}.Corr.C{2}=rec{ncg}.C;
       
        rec{ncg}.Corr.C{2}=multDimSum(rec{ncg}.Corr.C{2},4);          
        rec{ncg}.Corr.C{2}=ridgeDetection(rec{ncg}.Corr.C{2},1,32);
    
        
        rec{ncg}.Corr.C{2}=repmat(rec{ncg}.Corr.C{2},[1 size(rec{ncg}.z,2)/size(rec{ncg}.C,2)]);
        rec{ncg}.z=bsxfun(@times,rec{ncg}.z,rec{ncg}.Corr.C{2});
    end    
    
    %NYQUIST GHOST ESTIMATION AND CORRECTION
    dims2Sort=find(ismember(rec{ncg}.Corr.z{1},rec{ncg}.Plan.NDims+(1:2)));
    if isfield(rec{ncg},'P') && isfield(rec{ncg}.Corr,'P') && size(rec{ncg}.P,rec{ncg}.Corr.P{1}(2))>1   
        if any(typ2Rec==3)
            %USING AN EXTERNAL FILE FOR CALIBRATION DATA
            if ~isempty(rec{ncg}.Par.Mine.diInfo) && size(rec{ncg}.Par.Mine.diInfo,2)>=6
                rec{ncg}.P=dynInd(rec{ncg}.z,rec{ncg}.Par.Mine.diInfo(voInds,6)==0,11);
                if size(rec{ncg}.Par.Mine.diInfo,2)==7%Dirt readout flip due to inconsistencies in calibration in the DWI patch
                    for p=1:size(rec{ncg}.P,11)
                        if rec{ncg}.Par.Mine.diInfo(p,7);rec{ncg}.P=dynInd(rec{ncg}.P,p,11,conj(flip(dynInd(rec{ncg}.P,p,11),1)));end
                    end
                end
                rec{ncg}.z=dynInd(rec{ncg}.z,rec{ncg}.Par.Mine.diInfo(voInds,6)==1,11);
                voInds=voInds(rec{ncg}.Par.Mine.diInfo(voInds,6)==1);
                %for p=1:NPE;indPEV{p}=find(rec{ncg}.Par.Mine.diInfo(voInds,5)==rec{ncg}.Par.Mine.pedsUn(p));end
                for p=1:NPE;indPEV{p}=find(rec{ncg}.Par.Mine.diInfo(voInds,5)==p);end
            end
            %CHECK FOR PREVIOUS GHOSTING INFORMATION
            if size(rec{ncg}.Par.Mine.diInfo,2)<6 && rec{ncg}.Par.Mine.AdHocArray(1)==101 && rec{ncg}.Par.Mine.AdHocArray(4)>1 && rec{ncg}.Alg.GhosCorRef~=0 && rec{ncg}.Alg.UsePrevGhos %&& ~(isfield(rec{ncg}.Par.Mine,'StrFactorMax') && ~isempty(dims2Sort))
                indREF=find(rec{ncg}.Names.prot.A_Modal(1:rec{ncg}.Names.ind-1)==rec{ncg}.Par.Mine.Modal);   
                if ~isempty(indREF)%Candidate references
                    indPrev=length(indREF);                    
                    while indPrev>0%Check that is a valid reference
                        targetFile=strtrim(rec{ncg}.Names.prot.B_FileName(indREF(indPrev),:));                                                   
                        fileGHOST=fullfile(rec{ncg}.Names.pathOu,numbe2Modal(rec{ncg}.Par.Mine.Modal),sprintf('%s_Ny%s.mat',targetFile,rec{ncg}.Plan.Suff));
                        if exist(fileGHOST,'file')
                            warning('off','MATLAB:load:variableNotFound')
                            load(fileGHOST,'P','PExtr');
                            warning('on');
                            %if length(rec{ncg}.Par.Mine.pedsUn)==1 && size(P,11)>1;P=dynInd(P,rec{ncg}.Par.Mine.pedsUn,11);end
                            rec{ncg}.Corr.P{2}=P;
                            if exist('PExtr','var');rec{ncg}.Corr.PExtr{2}=PExtr;end
                            rec{ncg}.Alg.GhosCorRef=3;%Delays correction for ghosting. Uses SB data for correction
                            break
                        end
                        indPrev=indPrev-1;
                    end
                end
            end
            
            if rec{ncg}.Alg.GhosCorRef~=3
                %POLARITIES
                indE=cell(1,length(rec{ncg}.Plan.Par2Rec{7}));
                if rec{ncg}.Alg.GhosCorRef~=3;rec{ncg}.Corr.P{2}=rec{ncg}.P;end
                for e=rec{ncg}.Plan.Par2Rec{7}'%Different echoes
                    for n=typ2Rec';datTyp=rec{ncg}.Plan.Types{n};
                        if numel(rec{ncg}.Assign.z{7})>1;indE{e+1}=find(rec{ncg}.Assign.z{7}==e)';DimEcho=rec{ncg}.Corr.(datTyp){1}(7);else indE{e+1}=e+1;DimEcho=7;end
                        if ~ismember(n,4:5)
                            if max(indE{e+1})>size(rec{ncg}.Corr.(datTyp){5},DimEcho);fprintf('Problems in sorting, code needs extension to reconstruct scan %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end%Case 2014_09_10/SU_18303/su_10092014_1427347_15_2_dev1t2sdelayedrecsenseV4.raw
                            aux=dynInd(rec{ncg}.Corr.(datTyp){5},indE{e+1},DimEcho);
                            auxb.(datTyp){e+1}=gather(mapMat(aux,rec{ncg}.Corr.(datTyp){1}(2)));                        
                            auxb.(datTyp){e+1}=unique(auxb.(datTyp){e+1},'rows');
                            assert(size(auxb.(datTyp){e+1},1)<=2,'Polarity information of calibration not consistent among different EPI samples in typ %s for file %s',datTyp,name);aux=[];  
                        end
                    end                  
                    if any(abs(diff(auxb.P{e+1}))~=2);fprintf('The polarities are not alternating for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
                    if length(auxb.z{e+1})==length(auxb.P{e+1})
                        assert(gather(all(auxb.z{e+1}==auxb.P{e+1})),'The polarities are not the same in calibration and data for file %s',name);
                    end

                    %NYQUIST GHOST ESTIMATION
                    vshiftPlus=zeros(1,ndims(rec{ncg}.P));vshiftPlus(rec{ncg}.Corr.P{1}(2))=1;
                    vshiftMinu=vshiftPlus;vshiftMinu(rec{ncg}.Corr.P{1}(2))=-1;
                    if length(indE{e+1})>1
                        Paux=dynInd(rec{ncg}.P,indE{e+1},rec{ncg}.Corr.P{1}(7));
                        rec{ncg}.Corr.P{2}=dynInd(rec{ncg}.Corr.P{2},indE{e+1},rec{ncg}.Corr.P{1}(7),(Paux.^2.*conj(circshift(Paux,vshiftPlus).*circshift(Paux,vshiftMinu))).^(1/4));
                        indE{e+1}([1 end])=[];
                    else     
                        auxPa=(rec{ncg}.P.*conj(circshift(rec{ncg}.P,vshiftPlus))).^(1/2);
                        auxPb=(rec{ncg}.P.*conj(circshift(rec{ncg}.P,vshiftMinu))).^(1/2);
                        rec{ncg}.Corr.P{2}=(auxPa.*auxPb).^(1/2);
                        auxPa=[];auxPb=[];
                        %This operation may provoke numerical issues
                        %rec{ncg}.Corr.P{2}=(rec{ncg}.P.^2.*conj(circshift(rec{ncg}.P,vshiftPlus).*circshift(rec{ncg}.P,vshiftMinu))).^(1/4);
                    end      
                end
                rec{ncg}.Corr.P{2}=mean(rec{ncg}.Corr.P{2},4); 
                if numel(rec{ncg}.Assign.z{7})>1;indN=cat(2,indE{:});else indN=rec{ncg}.Assign.z{rec{ncg}.Corr.z{1}(2)};end
                for m=[2 5]
                    if rec{ncg}.Corr.z{1}(2)==2;indN=2:size(rec{ncg}.Corr.P{m},rec{ncg}.Corr.P{1}(2))-1;end     
                    rec{ncg}.Corr.P{m}=dynInd(rec{ncg}.Corr.P{m},indN,rec{ncg}.Corr.P{1}(2));     
                end

                %WE AVERAGE OVER LINES AND ECHOES
                if numel(rec{ncg}.Assign.z{7})>1
                    for m=1:length(nDirs);n=nDirs(m);
                        x{m}=dynInd(rec{ncg}.Corr.P{2},rec{ncg}.Corr.P{5}==n,rec{ncg}.Corr.P{1}(2));
                        x{m}=multDimMea(x{m},[rec{ncg}.Corr.P{1}(2) 7]);
                    end     
                else
                    for m=1:length(nDirs);n=nDirs(m);
                        for e=1:size(rec{ncg}.Corr.P{5},7)
                            y{e}=dynInd(rec{ncg}.Corr.P{2},{dynInd(rec{ncg}.Corr.P{5},e,7)==n,e},[rec{ncg}.Corr.P{1}(2) 7]);
                            y{e}=mean(y{e},rec{ncg}.Corr.P{1}(2));
                        end
                        x{m}=cat(7,y{:});
                        x{m}=mean(x{m},7);
                    end            
                end
                rec{ncg}.Corr.P{2}=horzcat(x{1},x{2});x=[];
                
                %WE EXTRAPOLATE IN THE SLICE DIRECTION---IT MAY BE BETTER TO USE PAPOULIS EXTRAPOLATION
                NPF=size(rec{ncg}.Corr.P{2});
                NPFR=NPF;NPFR(8)=2*NPF(8);
                rec{ncg}.Corr.PExtr{2}=resampling(rec{ncg}.Corr.P{2},NPFR,2);
                
                %WE FILTER IN THE SLICE DIRECTION---WE ALWAYS FILTER THE EXTRAPOLATED INFORMATION                  
                cosD=1;
                if rec{ncg}.Alg.GhosCorFilter                   
                    onePF=ones(1,max(numDims(rec{ncg}.Corr.P{2}),8));onePF(8)=NPF(8);                                                       
                    HPF=buildFilter((1+cosD)*onePF,'tukey',0.125,gpu,1,cosD);
                    rec{ncg}.Corr.P{2}=filtering(rec{ncg}.Corr.P{2},HPF,cosD);
                end
                onePF=ones(1,max(numDims(rec{ncg}.Corr.PExtr{2}),8));onePF(8)=2*NPF(8);
                HPF=buildFilter((1+cosD)*onePF,'tukey',0.25,gpu,1,cosD);
                rec{ncg}.Corr.PExtr{2}=filtering(rec{ncg}.Corr.PExtr{2},HPF,cosD);  
            end
        end
        
        %NYQUIST GHOST CORRECTION
        if ismember(rec{ncg}.Alg.GhosCorRef,[1 2])
            if rec{ncg}.Alg.GhosCorRef==1
                rec{ncg}.Corr.P{2}=ridgeDetection(rec{ncg}.Corr.P{2},1,32);
                rec{ncg}.Corr.PExtr{2}=ridgeDetection(rec{ncg}.Corr.PExtr{2},1,32);
            end  
            if rec{ncg}.Alg.parV.visGhost;visGhost(rec{ncg}.Corr.P{2});end
            if NPE~=1 && size(rec{ncg}.Par.Mine.diInfo,2)<6 && size(rec{ncg}.Corr.P{2},11)==1;rec{ncg}.Corr.P{2}=repmat(rec{ncg}.Corr.P{2},[ones(1,10) NPE]);end
            for m=1:length(nDirs);n=nDirs(m);
                for e=1:size(rec{ncg}.z,7)
                    indPE=dynInd(rec{ncg}.Corr.z{5},e,7)==n;
                    for p=1:size(rec{ncg}.Corr.P{2},11)
                        %x=dynInd(rec{ncg}.z,{indPE,e,indPEV{p}},[rec{ncg}.Corr.z{1}(2) 7 11]);                        
                        x=dynInd(rec{ncg}.z,{indPE,e,indPEV{rec{ncg}.Par.Mine.pedsUn(p)}},[rec{ncg}.Corr.z{1}(2) 7 11]);
                        if size(x,8)~=size(rec{ncg}.Corr.P{2},8);fprintf('Size of data and calibration not matched, the scan may have been interrupted for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end 
                        x=bsxfun(@times,x,exp(-1i*angle(dynInd(rec{ncg}.Corr.P{2},[m p],[2 11]))));
                        %rec{ncg}.z=dynInd(rec{ncg}.z,{indPE,e,indPEV{p}},[rec{ncg}.Corr.z{1}(2) 7 11],x);
                        rec{ncg}.z=dynInd(rec{ncg}.z,{indPE,e,indPEV{rec{ncg}.Par.Mine.pedsUn(p)}},[rec{ncg}.Corr.z{1}(2) 7 11],x);
                    end
                end
            end
        end
    end

    %ASSIGNMENT OF RECONSTRUCTION GRID
    %rec{ncg}=assignReconGrid(rec{ncg},dims2Sort,Enc);%This would replace the block of code below
    %if rec{ncg}.Fail==1;recOu=rec{ncg};return;end    
    
    N=size(rec{ncg}.z);    
    if isempty(dims2Sort);N(2:3)=Enc.AcqSize(2:3);
    elseif any(dims2Sort==12);N(2)=prod(Enc.AcqSize(2:3))*length(unique(rec{ncg}.Assign.z{12}(:)));
    else N(2)=prod(Enc.AcqSize(2:3));
    end
    N(end+1:6)=1;N(7)=rec{ncg}.Par.Encoding.NrEchoes;N(min(end,rec{ncg}.Plan.NDims)+1:rec{ncg}.Plan.NDims+2)=1;
    rec{ncg}.y=zeros(N,'like',rec{ncg}.z);

    if ~isempty(dims2Sort)    
        if isequal(dims2Sort,[2 7])
            for n=1:N(7)            
                assPE=rec{ncg}.Assign.z{2}(rec{ncg}.Assign.z{7}(:)==n-1);
                assert(numel(assPE)==diff(Enc.kRange{2}(n,:))+1,'The number of samples (%d) does not correspond to the prescribed ranges (%d) in echo %d for file %s',numel(assPE),diff(Enc.kRange{2}(n,:))+1,n,name);
                assert(all(ismember(assPE,Enc.kRange{2}(n,1):Enc.kRange{2}(n,2))),'Some data points are sampled outside the prescribed spectral region for file %s',name);
                assert(numel(unique(assPE))==numel(assPE),'Some data points have been sampled more than once or labeling may be corrupted for file %s',name);
                rec{ncg}.y=dynInd(rec{ncg}.y,{assPE(:)-min(Enc.kRange{2}(:,1))+1,n},[2 7],resPop(dynInd(rec{ncg}.z,rec{ncg}.Assign.z{7}(:)==n-1,rec{ncg}.Corr.z{1}(2)),rec{ncg}.Corr.z{1}(2),[],2));
            end            
        elseif isequal(dims2Sort,[2 3]) || isequal(dims2Sort,[2 3 12])
            for n=2:3
                if ~all(ismember(rec{ncg}.Assign.z{n}(:),Enc.kRange{n}(1,1):Enc.kRange{n}(1,2)))
                    fprintf('Some data points are sampled outside the prescribed spectral region for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;
                end
            end            
            if isequal(dims2Sort,[2 3])
                indPE=rec{ncg}.Assign.z{2}(:)-Enc.kRange{2}(1)+(rec{ncg}.Assign.z{3}(:)-Enc.kRange{3}(1))*Enc.AcqSize(2)+1;
            else
                indPE=rec{ncg}.Assign.z{2}(:)-Enc.kRange{2}(1)+(rec{ncg}.Assign.z{3}(:)-Enc.kRange{3}(1))*Enc.AcqSize(2)+rec{ncg}.Assign.z{12}(:)*prod(Enc.AcqSize(2:3))+1;
            end       
            if length(unique(indPE))~=length(indPE);fprintf('Some data points have been sampled more than once or labeling may be corrupted for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
            perm=1:rec{ncg}.Plan.NDims+2;perm(2)=rec{ncg}.Corr.z{1}(2);perm(rec{ncg}.Corr.z{1}(2))=2;
            if isequal(dims2Sort,[2 3])
                rec{ncg}.y=resPop(dynInd(rec{ncg}.y,indPE,2,permute(rec{ncg}.z,perm)),2,Enc.AcqSize(2:3),[2 3]);   
            else
                rec{ncg}.y=resPop(dynInd(rec{ncg}.y,indPE,2,permute(rec{ncg}.z,perm)),2,[Enc.AcqSize(2:3) N(2)/prod(Enc.AcqSize(2:3))],[2 3 12]); 
            end
            for d=dims2Sort             
                perm=1:rec{ncg}.Plan.NDims+2;perm(d)=rec{ncg}.Corr.z{1}(d);perm(rec{ncg}.Corr.z{1}(d))=d;             
                rec{ncg}.Assign.z{d}=permute(rec{ncg}.Assign.z{d},perm);
            end
            rec{ncg}.Corr.z{1}(dims2Sort)=dims2Sort;
        else
            fprintf('The conflict to be resolved (among dims%s) is not contemplated by the code\n',sprintf(' %d',dims2Sort));rec{ncg}.Fail=1;recOu=rec{ncg};return;
        end
    else        
        for n=1:N(7)
            for d=2:3
                assert(all(ismember(rec{ncg}.Assign.z{d}(:),Enc.kRange{d}(n,1):Enc.kRange{d}(n,2))),'Some data points are sampled outside the prescribed spectral region in echo %d and PE %d for file %s',n,d,name);
                if Enc.kRange{d}(n,2)-Enc.kGrid{d}(1)+1>size(rec{ncg}.z,d);fprintf('An interrupted sequence of data may have been encountered for dimension %d for file %s. The data will not be reconstructed.\n',d,name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
            end                                    
            rec{ncg}.y=dynInd(rec{ncg}.y,{(Enc.kRange{2}(n,1):Enc.kRange{2}(n,2))-Enc.kGrid{2}(1)+1,(Enc.kRange{3}(n,1):Enc.kRange{3}(n,2))-Enc.kGrid{3}(1)+1,n},[2 3 7],dynInd(rec{ncg}.z,n,7));
        end    
    end 
    rec{ncg}=rmfield(rec{ncg},'z');
    rec{ncg}.Dyn.Typ2Rec(rec{ncg}.Dyn.Typ2Rec==1)=[];
    

    %WE CROP THE SAMPLES FOR EXTENDED PHASE ENCODING IN DUAL SPIN-GRADIENT ECHO
    %TODO: ONCE THE ACQUISITION PATCH LABELLING GETS RIGHT THIS SHOULD NO LONGER BE NECESSARY
    if isfield(rec{ncg}.Par.Mine,'StrFactorMax') && ~isempty(dims2Sort)
        N=size(rec{ncg}.y,2);
        shift=N-1-rec{ncg}.Par.Mine.StrFactorMax;
        RP=diff(Enc.kRange{2}(2,:))+1;
        Enc.kRange{2}(2,1)=Enc.kRange{2}(2,1)+shift;
        Enc.PartFourier(2,2)=diff(Enc.kRange{2}(2,:)+1)/RP;   
        rec{ncg}.y=dynInd(rec{ncg}.y,2,7,circshift(dynInd(rec{ncg}.y,2,7),shift,2));
        rec{ncg}.y=dynInd(rec{ncg}.y,{1:shift,1:2},[2 7],0);
        
    end
    %WE REVERSE THE SECOND ECHO FOR DOSE
    if isfield(rec{ncg}.Par.Mine,'Spin2Echo') && rec{ncg}.Par.Mine.Spin2Echo && Enc.rRange{2}(2,1,1)==Enc.rRange{2}(1,1,1);rec{ncg}.y=dynInd(rec{ncg}.y,2,7,flip(dynInd(rec{ncg}.y,2,7),2));end%Not sure if this is rightly considering PF in the second echo, the code for doing so is at preprocessingData of ReconstructionsDebug04, however this rec{ncg}s may not be needed any longer
    
    %BACKTRANSFORM THE PHASE ENCODES, SET THE SLICES IN THIRD DIMENSION, CREATE
    %CONTAINERS FOR ENCODING AND RECONSTRUCTED INFORMATION
    if rec{ncg}.Alg.PlugNoise==3;rec{ncg}.y=plugNoise(rec{ncg}.y);end
    for n=2:3        
        for p=1:NPE
            if ~isempty(indPEV{p})
                A=Enc.DFTMH{p}{n};
                %Here it is size(A,2) while before it was size(A,1) because of Partial Fourier being the most important factor to be discounted here 
                rec{ncg}.y=dynInd(rec{ncg}.y,{1:size(A,1),indPEV{p}},[n 11],aplGPU(A,dynInd(rec{ncg}.y,{1:size(A,2),indPEV{p}},[n 11]),n)/sqrt(trace(A*A')/size(A,2)));
            end
        end
        if isempty(rec{ncg}.y);fprintf('Data appears empty for file %s\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end%This error is observed in usual corruption of spin-echo 2 phase encode of dHCP, example in series 6/8 of 2019_09_18/dh_128031
        rec{ncg}.y=dynInd(rec{ncg}.y,1:size(A,1),n);    
    end    
    if rec{ncg}.Alg.PlugNoise==4;rec{ncg}.y=plugNoise(rec{ncg}.y);end    
    rec{ncg}.Dyn.Typ2Rec=vertcat(rec{ncg}.Dyn.Typ2Rec,6);
    N=size(rec{ncg}.y);N(end+1:rec{ncg}.Plan.NDims)=1;    
    
    %GIBBS RINGING FILTER FOR COIL INFORMATION
    if rec{ncg}.Par.Mine.Modal==2
        sizFilt=N(1:3);
        sp=Enc.AcqVoxelSize;
        sp(3)=sp(3)+rec{ncg}.Par.Labels.SliceGaps(1);
        spaFilt=max(sp)./sp;
        spaFilt=spaFilt.*[1 max(rec{end}.Assign.z{2})/max(Enc.kRange{2}(:,2)) max(rec{end}.Assign.z{3})/max(Enc.kRange{3}(:,2))];
        spaFilt(2:3)=spaFilt(2:3)*rec{ncg}.Alg.parS.ResolRatio(ncg);
        %spaFilt(2:3)=rec{ncg}.Alg.parS.ResolRatio(ncg);
        H=buildFilter(2*sizFilt,'tukeyIso',spaFilt,gpu,rec{ncg}.Alg.parS.GibbsRingi(ncg),1);         
        rec{ncg}.y=filtering(rec{ncg}.y,H,1);    
    end
    
    if all(N([3 8])~=1);fprintf('Both the second phase encode and the slice direction are non-singletons for file %s. A reconstruction method has to be put in place for this\n',name);rec{ncg}.Fail=1;recOu=rec{ncg};return;end
    if N(8)>1 && N(3)==1;perm=1:rec{ncg}.Plan.NDims;perm([3 8])=[8 3];rec{ncg}.y=permute(rec{ncg}.y,perm);end
    rec{ncg}.Enc=Enc;    
    for n=6;datTyp=rec{ncg}.Plan.Types{n};
        if rec{ncg}.Dyn.Debug>=2;printResults;end
    end    
end
%COMBINE COIL RESULTS
recOu=rec{1};
if NCG==2
    recOu.Par.Labels.CoilNrsPerStack{2}=rec{2}.Par.Labels.CoilNrsPerStack{1};%As for the body coil we assume the order does not matter
    recOu.y=cat(4,recOu.y,rec{2}.y);
    if length(unique(vertcat(recOu.Par.Labels.CoilNrsPerStack{:})))==length(vertcat(recOu.Par.Labels.CoilNrsPerStack{:}))
        recOu.y=dynInd(recOu.y,vertcat(recOu.Par.Labels.CoilNrsPerStack{:})'+1,4,recOu.y);
    end
end

%TEST OF COMMON ENCODINGS BETWEEN ECHOES
if size(unique(recOu.Enc.kRange{2},'rows'),1)~=1;fprintf('The encoding is different for different echoes:%s, note this is only supported by the SB reconstructor\n',sprintf(' %d',recOu.Enc.kRange{2}));end

%MULTIBAND ENCODING
AdHocArray=recOu.Par.Mine.AdHocArray;
if isempty(AdHocArray) || all(AdHocArray==0:127);recOu.Enc.MB=1;return;end
assert(AdHocArray(1)~=102,'SIR reconstruction is not supported any longer'); 
if AdHocArray(1)==101;recOu.Enc.MB=AdHocArray(4);else recOu.Enc.MB=1;end
if recOu.Enc.MB==1 && AdHocArray(7)<=1 && ((recOu.Par.Mine.Modal==9 && ~strcmp(recOu.Par.Scan.Technique,'SEEPI')) || recOu.Par.Mine.Modal==10)
    if isfield(recOu.Corr,'P') && size(recOu.Corr.P{2},2)==2;recOu.Dyn.Typ2Wri(3)=1;end
    return;
end

%Z-BLIP ENCODING STRENGTHS
if isempty(recOu.Par.Encoding.KzRange)
    recOu.Enc.RecSize(3)=size(recOu.y,3);
    recOu.Enc.FOVSize(3)=size(recOu.y,3)*recOu.Enc.MB;
    recOu.Enc.SlicDist=recOu.Par.Labels.VoxelSizes(3)+recOu.Par.Labels.SliceGaps(1);%Slice separation
    if AdHocArray(3)==3%TO DO---USE A GENERATEGRID FOR THIS
        sliceGap=AdHocArray(2)/recOu.Enc.SlicDist;
        recOu.Enc.ExcGrid=(0:recOu.Enc.RecSize(3)-1)';
        ExcGridOffs=(0:recOu.Enc.MB-1)*sliceGap;
        ExcGridOffs=ExcGridOffs-(recOu.Enc.MB-1)*(sliceGap-recOu.Enc.RecSize(3))/2;
        recOu.Enc.ExcGrid=bsxfun(@plus,recOu.Enc.ExcGrid,ExcGridOffs);
        recOu.Enc.ExcGrid=reshape(recOu.Enc.ExcGrid,[1 recOu.Enc.FOVSize(3)]);
    else
        recOu.Enc.ExcGrid=0:recOu.Enc.FOVSize(3)-1;
    end
end

%TEST OF COMMON ENCODINGS BETWEEN ECHOES
if size(unique(recOu.Enc.kRange{2},'rows'),1)~=1;fprintf('The encoding is different for different echoes:%s, not supported by the MB reconstructor\n',sprintf(' %d',recOu.Enc.kRange{2}));recOu.Fail=1;return;end

%FREQUENCY STABILIZATION (activated ony if fMRI and fetal cases)
if ( ~isempty(regexp(recOu.Names.Name,'sf','ONCE')) || ~isempty(regexp(recOu.Names.Name,'fd','ONCE')) || ~isempty(regexp(recOu.Names.Name,'f[5-9]','ONCE')) ) && recOu.Par.Mine.Modal==9 && recOu.Alg.StabilFreq
    if recOu.Dyn.Batch==1;recOu.yf=dynInd(recOu.y,1,5);end
    perm=1:rec{ncg}.Plan.NDims;perm(4:5)=[5 4];
    for n=1:size(recOu.y,5)
        T=computeFrequencyStabilization(permute(cat(5,recOu.yf,dynInd(recOu.y,n,5)),perm),2,[],2,recOu);
        %fprintf('Stabilization shift: %.4f\n',T{2});
        [y,ESFOV,TE]=applyFrequencyStabilization(dynInd(recOu.y,n,5),T,recOu);
        recOu.y=dynInd(recOu.y,n,5,y);        
        F=convertB0Field(gather(T{2}),TE,ESFOV,'pix','Hz',size(recOu.y,2));
        if ~isfield(recOu,'F');recOu.F=F;else recOu.F=cat(4,recOu.F,F);end
        recOu.Dyn.Typ2Wri(23)=1;
    end
end

function printResults
    fprintf('Size of data typ %s:%s\n',datTyp,sprintf(' %d',size(rec{ncg}.(datTyp))));      
end

end
