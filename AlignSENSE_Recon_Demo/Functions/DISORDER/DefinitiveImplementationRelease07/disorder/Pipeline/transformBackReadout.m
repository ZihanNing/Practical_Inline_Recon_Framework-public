function rec=transformBackReadout(rec,Enc,indPEV,full)

%TRANSFORMBACKREADOUT   Performs the basic readout Fourier inversion of acquired data
%   REC=TRANSFORMBACKREADOUT(REC)
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), the
%   .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)) and the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types))
%   * FULL transforms to the full FOV in the readout direction
%   * REC is a reconstruction structure with backtransformed data (rec.(rec.PlanTypes)) and encoding (rec.Enc) information
%

if nargin<4 || isempty(full);full=0;end

typ2Rec=rec.Dyn.Typ2Rec;name=rec.Names.Name;gpu=rec.Dyn.GPU;
nDirs=unique(rec.Corr.z{5}(:))';
NPE=length(indPEV);

%ACTUAL TRANSFORM
if rec.Dyn.GPU==2;rec.N=gpuArray(rec.N);end%From saving memory
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if n~=5
        if rec.Dyn.GPU==2;rec.(datTyp)=gpuArray(rec.(datTyp));end%From saving memory
        dims=1:rec.Plan.NDims+2;
        if rec.Corr.(datTyp){1}(2)==0 || numel(rec.Corr.(datTyp){5})/size(rec.Corr.(datTyp){5},4)==1
            if rec.Dyn.Debug>=1;fprintf('Unclear interpretation of calibration as only a single PE is present in typ %s for file %s\n',datTyp,name);end
            rec.Corr.(datTyp){1}(2)=2;
        end

        %NOISE STANDARDIZATION
        if rec.Alg.NoiseStand && rec.Par.Mine.Modal~=2;rec.(datTyp)=standardizeCoils(rec.(datTyp),rec.N);end

        %TRANSFORM BACK THE READOUT
        aux=gather(mapMat(rec.Corr.(datTyp){5},rec.Corr.(datTyp){1}(2)));
        assert(size(unique(aux,'rows'),1)<=2,'Polarity information of calibration not consistent among different EPI samples in typ %s for file %s',datTyp,name);aux=[]; 

        dimsInt=setdiff(dims,rec.Corr.(datTyp){1}([2 7]));
        rec.Corr.(datTyp){5}=dynInd(rec.Corr.(datTyp){5},ones(1,length(dimsInt)),dimsInt);

        if ~(rec.Corr.(datTyp){1}(2)<7);E=0;else E=rec.Plan.Par2Rec{7}';end
        if n==4;E=0:size(rec.(datTyp),7)-1;end
        for e=E%Different echoes
            ex=ones(1,rec.Plan.NDims+1);ex(7-single((rec.Corr.(datTyp){1}(2))<7))=e+1;   
            aux=dynInd(rec.Corr.(datTyp){5},ex,setdiff(dims,rec.Corr.(datTyp){1}(2)));                             
            for nd=nDirs
                for p=1:NPE
                    %if n==1;indV=indPEV{p};elseif rec{ncg}.Par.Mine.pedsUn(p)==rec{ncg}.Par.Mine.Nat || length(rec{ncg}.Par.Mine.pedsUn)==1;indV=1;else indV=[];end
                    if ismember(n,[1 4]);indV=indPEV{p};elseif p==rec.Par.Mine.Nat;indV=1;else indV=[];end
                    if n==4;indV=indV(1:min(length(indV),size(rec.(datTyp),11)));end
                    if ~isempty(indV)
                        x=dynInd(rec.(datTyp),{aux(:)==nd,e+1,indV},[rec.Corr.(datTyp){1}(2) 7 11]);
                        if ~isempty(x)
                            N=size(x);
                            if ~full;N(1)=Enc.RecSize(1);else N(1)=Enc.AcqSize(1);end
                            if rec.Alg.PlugNoise==1 && n==1;x=plugNoise(x);end
                            if ~full
                                if nd>0;A=Enc.DFTMH{p}{1};else A=conj(Enc.DFTMH{p}{1});end     
                            else
                                if nd>0;A=Enc.DFTMHAcq{p}{1};else A=conj(Enc.DFTMHAcq{p}{1});end
                            end
                            x=aplGPU(A,x,1)/sqrt(trace(A*A')/size(A,1));
                            if rec.Alg.PlugNoise==2 && n==1;x=plugNoise(x);end
                            rec.(datTyp)=dynInd(rec.(datTyp),{1:N(1),aux(:)==nd,e+1,indV},[1 rec.Corr.(datTyp){1}(2) 7 11],x);
                        end
                    end
                end
            end
        end
        x=[];rec.(datTyp)=dynInd(rec.(datTyp),1:N(1),1);
    end
end
