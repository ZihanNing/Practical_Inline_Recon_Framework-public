function rec=assignReconGrid(rec,dims2Sort,Enc)

%ASSIGNRECONGRID   Assigns the reconstruction grid
%   REC=ASSIGNRECONGRID(REC)
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), the
%   .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)) and the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types))
%   * FULL transforms to the full FOV in the readout direction
%   * REC is a reconstruction structure with backtransformed data (rec.(rec.PlanTypes)) and encoding (rec.Enc) information
%

N=size(rec.z);    
if isempty(dims2Sort);N(2:3)=Enc.AcqSize(2:3);
elseif any(dims2Sort==12);N(2)=prod(Enc.AcqSize(2:3))*length(unique(rec.Assign.z{12}(:)));
else N(2)=prod(Enc.AcqSize(2:3));
end
N(end+1:6)=1;N(7)=rec.Par.Encoding.NrEchoes;N(min(end,rec.Plan.NDims)+1:rec.Plan.NDims+2)=1;
rec.y=zeros(N,'like',rec.z);

if ~isempty(dims2Sort)    
    if isequal(dims2Sort,[2 7])
        for n=1:N(7)            
            assPE=rec.Assign.z{2}(rec.Assign.z{7}(:)==n-1);
            assert(numel(assPE)==diff(Enc.kRange{2}(n,:))+1,'The number of samples (%d) does not correspond to the prescribed ranges (%d) in echo %d for file %s',numel(assPE),diff(Enc.kRange{2}(n,:))+1,n,name);
            assert(all(ismember(assPE,Enc.kRange{2}(n,1):Enc.kRange{2}(n,2))),'Some data points are sampled outside the prescribed spectral region for file %s',name);
            assert(numel(unique(assPE))==numel(assPE),'Some data points have been sampled more than once or labeling may be corrupted for file %s',name);
            rec.y=dynInd(rec.y,{assPE(:)-min(Enc.kRange{2}(:,1))+1,n},[2 7],resPop(dynInd(rec.z,rec.Assign.z{7}(:)==n-1,rec.Corr.z{1}(2)),rec.Corr.z{1}(2),[],2));
        end            
    elseif isequal(dims2Sort,[2 3]) || isequal(dims2Sort,[2 3 12])
        for n=2:3
            if ~all(ismember(rec.Assign.z{n}(:),Enc.kRange{n}(1):Enc.kRange{n}(2)))
                fprintf('Some data points are sampled outside the prescribed spectral region for file %s\n',name);rec.Fail=1;return
            end
        end
        if isequal(dims2Sort,[2 3])
            indPE=rec.Assign.z{2}(:)-Enc.kRange{2}(1)+(rec.Assign.z{3}(:)-Enc.kRange{3}(1))*Enc.AcqSize(2)+1;
        else
            indPE=rec.Assign.z{2}(:)-Enc.kRange{2}(1)+(rec.Assign.z{3}(:)-Enc.kRange{3}(1))*Enc.AcqSize(2)+rec.Assign.z{12}(:)*prod(Enc.AcqSize(2:3))+1;
        end        
        if length(unique(indPE))~=length(indPE);fprintf('Some data points have been sampled more than once or labeling may be corrupted for file %s\n',name);rec.Fail=1;return;end
        perm=1:rec.Plan.NDims+2;perm(2)=rec.Corr.z{1}(2);perm(rec.Corr.z{1}(2))=2;
        if isequal(dims2Sort,[2 3])
            rec.y=resPop(dynInd(rec.y,indPE,2,permute(rec.z,perm)),2,Enc.AcqSize(2:3),[2 3]);   
        else
            rec.y=resPop(dynInd(rec.y,indPE,2,permute(rec.z,perm)),2,[Enc.AcqSize(2:3) N(2)/prod(Enc.AcqSize(2:3))],[2 3 12]); 
        end
        for d=dims2Sort             
            perm=1:rec.Plan.NDims+2;perm(d)=rec.Corr.z{1}(d);perm(rec.Corr.z{1}(d))=d;             
            rec.Assign.z{d}=permute(rec.Assign.z{d},perm);
        end
        rec.Corr.z{1}(dims2Sort)=dims2Sort;
    else
        fprintf('The conflict to be resolved (among dims%s) is not contemplated by the code\n',sprintf(' %d',dims2Sort));rec.Fail=1;return;
    end
else        
    for n=1:N(7)
        for d=2:3
            assert(all(ismember(rec.Assign.z{d}(:),Enc.kRange{d}(n,1):Enc.kRange{d}(n,2))),'Some data points are sampled outside the prescribed spectral region in echo %d and PE %d for file %s',n,d,name);
            if Enc.kRange{d}(n,2)-Enc.kGrid{d}(1)+1>size(rec.z,d);fprintf('An interrupted sequence of data may have been encountered for dimension %d for file %s. The data will not be reconstructed.\n',d,name);rec.Fail=1;return;end
        end                                    
        rec.y=dynInd(rec.y,{(Enc.kRange{2}(n,1):Enc.kRange{2}(n,2))-Enc.kGrid{2}(1)+1,(Enc.kRange{3}(n,1):Enc.kRange{3}(n,2))-Enc.kGrid{3}(1)+1,n},[2 3 7],dynInd(rec.z,n,7));
    end    
end 
rec=rmfield(rec,'z');
rec.Dyn.Typ2Rec(rec.Dyn.Typ2Rec==1)=[];