function rec=writeData(rec,traMR)

%WRITEDATA   Writes reconstructed data to disk and assigns geometry
%information
%   REC=WRITEDATA(REC,{TRAMR})
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   * {TRAMR} indicates wheter to perform a transformation from MPS to REC 
%   coordinates, defaults to 0
%   ** REC is a reconstruction structure with reconstructed data and
%   surrogates information (rec.(rec.Plan.Types)) and removed information 
%   for correction (rec.Corr.(rec.Plan.Types)) and sorting 
%   (rec.Assign(rec.Plan.Types))
%

if nargin<2 || isempty(traMR);traMR=0;end

%GENERATE THE FILE NAME
fileName=generateNIIFileName(rec);

%WRITE JSON
if rec.Alg.OnlyJSON;writeJSON(rec,fileName);end

if rec.Fail || ~any(rec.Dyn.Typ2Wri)
    if nargout==0;rec=[];end%To save memory
    return;
end
rec.Dyn.Typ2Rec(ismember(rec.Dyn.Typ2Rec,[3 6]))=[];%YB: so 3 and 6 (EPIread and spectra) not to write

typ2Rec=rec.Dyn.Typ2Rec;

% %BACK TO ORIGINAL ROI
% if isfield(rec.Enc,'ROI')
%     for n=typ2Rec';datTyp=rec.Plan.Types{n};
%         if rec.Dyn.Typ2Wri(n) && ~ismember(n,24:25);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0);end
%     end
% end
%rec.Enc.ROI
%pause

%ROTATE THE RECONSTRUCTIONS TO SYSTEM INDEPENDENT OF ENCODING AXES
if traMR;rec=mps2Rec(rec);end

if ~rec.Par.Mine.Proce
    %REMOVE THE OVERDECODING FROM THOSE DATASETS THAT NEED TO
    for n=typ2Rec';datTyp=rec.Plan.Types{n};    
        if ~ismember(n,[5 9 10 12 16]);rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);end
    end
    if rec.Alg.EstimGFact(1) && isfield(rec,'G') && isfield(rec,'M') && rec.Dyn.Typ2Wri(9)
        if ~ismember(7,typ2Rec');rec.S=removeOverencoding(rec.S,rec.Alg.OverDec);end
        if ~ismember(27,typ2Rec') && isfield(rec,'W');rec.W=removeOverencoding(rec.W,rec.Alg.OverDec);end
        if ~ismember(8,typ2Rec');rec.M=removeOverencoding(rec.M,rec.Alg.OverDec);end
        perm=1:6;perm(4:6)=[5 4 6];
        if size(rec.M,5)==4;M=permute(dynInd(rec.M,rec.Par.Mine.pedsUn,5),perm);
        elseif size(rec.M,5)==length(rec.Par.Mine.pedsUn);M=permute(rec.M,perm);
        else error('Unknown PE configuration at writeData. Number of reformatted masks: %d / Number of PEs: %d',size(rec.M,5),length(rec.Par.Mine.pedsUn));
        end  
        if size(rec.S,5)==4;S=permute(dynInd(rec.S,rec.Par.Mine.pedsUn,5),perm);
        elseif size(rec.S,5)==length(rec.Par.Mine.pedsUn);S=permute(rec.S,perm);
        else error('Unknown PE configuration at writeData. Number of reformatted masks: %d / Number of PEs: %d',size(rec.S,5),length(rec.Par.Mine.pedsUn));
        end  
        
        if rec.Alg.UseSoftMasking>0;M(:)=1;end
        S=dynInd(S,1,6);
        %%%%THIS LINE GIVES PROBLEMS IN CASE OF SEVERAL PES AND SEPARATE NOISE ESTIMATION    
        rec.G=cat(5,sqrt(rec.Enc.Accel)*rec.G,bsxfun(@times,bsxfun(@times,M,rec.G),sqrt(1./(sum(abs(S).^2,5)+1e-9))));
    end
    
    %FLIP GRIDS FOR MULTIPLE PHASE ENCODE DATASETS
    rec=PEregrid(rec);
    if rec.Fail
        if nargout==0;rec=[];end%To save memory
        return;
    end
end

%GENERATE THE GEOMETRICAL INFORMATION
if isfield(rec,'MT')%Case procePipelineTransform
    MT=rec.MT{1};
    MS=rec.MS{1};
else
    [MS,MT]=generateNIIInformation(rec.Par);
end
if isfield(rec,'Enc') && isfield(rec.Enc,'ROI')
    vROI=rec.Enc.ROI(:,1)-1;vROI(4)=1; 
    vROI=MT*vROI;
    MTROI=MT;
    MTROI(1:3,4)=vROI(1:3);
else
    vROI=zeros(4,1);vROI(4)=1;
    overDec=rec.Alg.OverDec;
    overDec(overDec>1)=1;overDec=abs(overDec);        
    if ~isempty(rec.Par.Mine.Signs);overDec(1:2)=max(overDec(1:2));end

    for n=1:3
        if overDec(n)>1
            N=size(rec.x,n);Nor=round(N/overDec(n));
            vROI(n)=vROI(n)-ceil((N-Nor)/2);
        end
    end
    vROI=MT*vROI;
    MT(1:3,4)=vROI(1:3);
end

%GENERATE STRUCTURES FOR WRITING AND WRITE IMAGES
c=1;
for n=typ2Rec';datTyp=rec.Plan.Types{n};datName=rec.Plan.TypeNames{n};    
    if rec.Dyn.Typ2Wri(n)
        if n==8;rec.(datTyp)=single(abs(rec.(datTyp)));end%Soft mask
        if n==7 && rec.Alg.NoiseStand && rec.Par.Mine.Modal~=2 && isfield(rec,'N');rec.(datTyp)=standardizeCoils(rec.(datTyp),rec.N,1);end%Coil de-standardization        
        if n==29 %YB added
            x{c}=dynInd( rec.(datTyp),1,6);MSV{c}=MS;
            x{c+1}=dynInd( rec.(datTyp),2,6);MSV{c+1}=MS;
            if isfield(rec,'Enc') && isfield(rec.Enc,'ROI');MTV{c}=MTROI;MTV{c+1}=MTROI;else MTV{c}=MT;MTV{c+1}=MT;end
        else 
            x{c}=rec.(datTyp);MSV{c}=MS;
        end
        
        if ~ismember(n,24:25) && n~= 29 && isfield(rec,'Enc') && isfield(rec.Enc,'ROI');MTV{c}=MTROI;else MTV{c}=MT;end
        if n==11 && rec.Par.Mine.Modal==4;datName='B1';end
        if n==29 %YB added
            suff{c}=strcat([datName '_1'],rec.Plan.Suff,rec.Plan.SuffOu);      
            suff{c+1}=strcat([datName '_2'],rec.Plan.Suff,rec.Plan.SuffOu);  
            c=c+2;
        else
            suff{c}=strcat(datName,rec.Plan.Suff,rec.Plan.SuffOu);  
            c = c+1;
        end
    end
end

if c>1;writeNII(fileName,suff,x,MSV,MTV);end
if rec.Dyn.Typ2Wri(12);writeJSON(rec,fileName);writeJSON(rec,fileName,1);end

%WRITE SURROGATE INFORMATION
if rec.Dyn.Typ2Wri(3)
        P=gather(rec.Corr.P{2});
    if isfield(rec.Corr,'PExtr')
        PExtr=gather(rec.Corr.PExtr{2});
        save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{3},rec.Plan.Suff,rec.Plan.SuffOu)),'P','PExtr');
    else
        save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{3},rec.Plan.Suff,rec.Plan.SuffOu)),'P');
    end
end
if rec.Dyn.Typ2Wri(17)
    T=gather(rec.T);
    if isfield(rec.Par.Mine,'DISORDER');MotionInfo=rec.Par.Mine.DISORDER;else;MotionInfo=[];end
    MotionInfo.timeState = gather(rec.timeState);
    MotionInfo.Par = rec.Par;
    MotionInfo.T = T;
    if isfield(rec,'E') && isfield(rec.E,'Dl') && isfield(rec.E.Dl,'Tr');MotionInfo.T_RAS = rec.E.Dl.Tr;end
    save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{17},rec.Plan.Suff,rec.Plan.SuffOu)),'T','MotionInfo');
end
%Save PT information
if isfield(rec,'PT')
    PT = gatherStruct(rec.PT);
    PT.pTimeRes=[];PT.pSliceRes=[];PT.pSliceImage=[];PT.pSlice=[];%PT.pTimeTest=[];%Can be found in raw data .mat file
    if rec.Dyn.Typ2Wri(17) && isfield( rec, 'T'); T = gather(rec.T);else; T=[];end
    if isfield(rec.Alg,'parXT') && isfield(rec.Alg.parXT,'PT'); PT.Info =  gatherStruct (rec.Alg.parXT.PT);end
    save(sprintf('%s_%s.mat',fileName,strcat('PTInfo',rec.Plan.Suff,rec.Plan.SuffOu)),'PT','T');
end
if isfield(rec,'E') && isfield(rec.E,'Db') %YB added to inspect basis coefficients later
    Db=gather(rec.E.Db);
    Db = rmfield(Db,'B');%avoid file to be big
    timeState = gather(rec.timeState);
    save(sprintf('%s_%s.mat',fileName,strcat('Db',rec.Plan.Suff,rec.Plan.SuffOu)),'Db','timeState');
end

if rec.Dyn.Typ2Wri(23)
    F=gather(rec.F);save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{23},rec.Plan.Suff,rec.Plan.SuffOu)),'F');
end

if rec.Dyn.Typ2Wri(25)
    t=cellFun(rec.t,@gather);save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{25},rec.Plan.Suff,rec.Plan.SuffOu)),'t');
end

if nargout==0;rec=[];end%To save memory