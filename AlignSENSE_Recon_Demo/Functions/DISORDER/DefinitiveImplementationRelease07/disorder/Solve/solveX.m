function rec=solveX(rec)

%SOLVEX   Estimates the reconstructed data
%   REC=SOLVEX(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated reconstructed data (rec.x)
%

gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU && rec.Dyn.GPU~=0);

Alg=rec.Alg;
parX=Alg.parX;
%NPE=size(rec.S,5);
NPE=length(rec.Par.Mine.pedsUn);
voInds=rec.Dyn.Ind2Rec{11}+1;%Volumetric indexes
if ~isempty(rec.Par.Mine.diInfo) && size(rec.Par.Mine.diInfo,2)>=6;voInds=voInds(rec.Par.Mine.diInfo(voInds,6)==1);end
gp=[];xp=[];cp=[];indPEV=cell(1,NPE);
sepno=(Alg.UnfoldFat && Alg.UseSoftMasking && Alg.EstimCFact) || strcmp(Alg.solverType,'IRWLS');%Flag to compute the noise volume by volume
unfat=(Alg.UnfoldFat && Alg.UseSoftMasking && Alg.EstimCFact && isfield(rec,'W'));

for p=1:NPE
    estGF=Alg.EstimGFact;        
    %Data
    if rec.Par.Mine.Modal==2;y=dynInd(rec.y,rec.Par.Labels.CoilNrsPerStack{1}+1,4);else y=rec.y;end
    %PEs of interest 
    if ~isempty(rec.Par.Mine.diInfo);indPEV{p}=find(rec.Par.Mine.diInfo(voInds,5)==rec.Par.Mine.pedsUn(p));else indPEV{p}=1:length(voInds);end    
    if ~isempty(indPEV{p})
        %[E,EH]=MBEncoding(rec,p);
        [E,EH,rec.Fail]=MBEncoding(rec,rec.Par.Mine.pedsUn(p));%This may need similar corrections as below for S and M in case it is applied for non-conventional protocols
        if rec.Fail==1;return;end
        if isfield(rec.Par.Mine,'ShiftIsoc') && Alg.JointSplitScans
            indSplit=rec.Par.Mine.splitInfo(voInds);
            indSplit=indSplit(indPEV{p});
            P.PV=E.Gh.Pf;P.PVH=EH.Gh.Pb;
            B.BV=E.Bf;B.BVH=EH.Bb;
        end
        if ~isfield(rec,'S');fprintf('Sensitivitites could not be identified\n');rec.Fail=1;return;end
        if size(rec.S,5)==4
            E.Sf=dynInd(rec.S,rec.Par.Mine.pedsUn(p),5);
            if isfield(rec,'W');E.WW=dynInd(rec.W,rec.Par.Mine.pedsUn(p),5);end
        elseif size(rec.S,5)==NPE
            E.Sf=dynInd(rec.S,p,5);
            if isfield(rec,'W');E.WW=dynInd(rec.W,p,5);end
            %E.WW=dynInd(E.WW,2:size(E.WW,6),6,0);
        else fprintf('Unknown PE configuration at solveX. Number of reformatted coil maps: %d / Number of PEs: %d\n',size(rec.S,5),NPE);rec.Fail=1;return;
        end        
        y=dynInd(y,indPEV{p},11);

        %Encoder/Decoder
        %Coil compression
        if Alg.parV.visSensitivities;visSensitivities(E.Sf,Alg.parV.visSensitivities,rec.Par.Scan.AcqVoxelSize);end
        if size(E.Sf,6)==1
            if rec.Par.Mine.Modal~=7;[E.Sf,y]=compressCoils(E.Sf,parX.perc,y,0);else [E.Sf,y]=compressCoils(E.Sf,parX.perc,y);end
        end
        if gpuIn;y=gpuArray(y);end
        NS=size(E.Sf);NS(end+1:4)=1;
        if rec.Dyn.Debug>=2 && ~isempty(parX.perc) && parX.perc(1)<1;fprintf('Number of compressed coil elements at %2.1f%%: %d\n',parX.perc(1)*100,NS(4));end
        if prod(NS(1:3))>rec.Dyn.MaxMem(2);E.bS=[1 1];else E.bS=[1 NS(4)];end%Block sizes
        E.dS=[1 NS(4)];E.oS=[1 1];%End and start index of operation of block size
        E.Uf=cell(1,3);for n=1:length(E.Uf);E.Uf{n}.NY=rec.Enc.RecSize(n);E.Uf{n}.NX=rec.Enc.FOVSize(n);end
        E.Je=0;%Joint encoding and decoding
        E.pe=Alg.PE;EH.pe=Alg.PE;%Phase encode dimension

        %GIBBS RINGING FILTERING
        if parX.UseGiRingi==1 && any(parX.GibbsRingi>0)    
            if rec.Enc.NDimsGR~=length(parX.GibbsRingi)      
                H=buildFilter(rec.Enc.RecSize(1:rec.Enc.NDimsGR),'tukeyIso',ones(1,rec.Enc.NDimsGR),gpuIn,parX.GibbsRingi(1));
            else     
                H=buildFilter(rec.Enc.RecSize(1:rec.Enc.NDimsGR),'tukey',ones(1,rec.Enc.NDimsGR),gpuIn,parX.GibbsRingi);
            end
            y=filtering(y,H);
        end
        
        %THIS DEFINES WHETHER THE RESIDUALS LIE IN FOURIER SPACE OR NOT
        E.Es=1;
        if (~isfield(E,'Es') || E.Es==0) && isfield(E,'Ef');y=aplGPU(E.Ef,y,2);end%For MB or ghosting correction, PE data should be in Fourier domain                
        if Alg.PlugNoise==5;y=plugNoise(y);end%To plug noise samples
        
        %Decoder
        EH.Ub=cell(1,3);for n=1:length(EH.Ub);EH.Ub{n}.NY=rec.Enc.RecSize(n);EH.Ub{n}.NX=rec.Enc.FOVSize(n);end
        [E.UAf,EH.UAb]=buildFoldM(rec.Enc.FOVSize(1:length(E.Uf)),rec.Enc.RecSize(1:length(E.Uf)),gpuIn);
        %Precondition
        if gpuIn && rec.Par.Mine.Modal==7;E.Sf=gpuArray(E.Sf);end
        A.Se=(sum(abs(conj(E.Sf).*E.Sf),4)+1e-9).^(-1);       
        %Constrain
        if size(rec.M,5)==4;C.Ma=dynInd(rec.M,rec.Par.Mine.pedsUn(p),5);
        elseif size(rec.M,5)==NPE;C.Ma=dynInd(rec.M,p,5);
        else fprintf('Unknown PE configuration at solveX. Number of reformatted masks: %d / Number of PEs: %d\n',size(rec.M,5),NPE);rec.Fail=1;return;
        end
        if gpuIn;A.Se=gpuArray(A.Se);C.Ma=gpuArray(C.Ma);end
        if Alg.UseSoftMasking==2;C.Ma(:)=1;end
        if E.Uf{3}.NX~=1 && size(C.Ma,3)~=E.Uf{3}.NX;fprintf('Acquisition may be incomplete as number of slices do not match for inverted (%d) and to be reconstructed (%d) FOVs\n',E.Uf{3}.NX,size(C.Ma,3));rec.Fail=1;return;end

        %Regularizer
        R=[];       
        if any(~ismember(C.Ma(:),[0 1]))
            %visReconstruction(C.Ma(:,:,end),0)
            M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpuIn,1,1);
            C.Ma=abs(filtering(C.Ma,M,1));       
            if size(E.Sf,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction--note this may introduce truncation artifacts, not seen so far
                NMa=[E.Uf{1}.NX E.Uf{2}.NX];
                for n=1:2
                    if ~isempty(E.Uf{n});C.Ma=fold(C.Ma,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
                end         
                C.Ma(:)=Alg.parE.lambda;
                C.Ma=resampling(C.Ma,NMa,2);
            end          
            if isfield(rec,'W')
                Waux=2*(E.WW-Alg.parE.eigSc(1))./(1-Alg.parE.eigSc(1))-1;
                Waux=0.5+0.5*tanh(3*Waux);%"Standard" sigmoid
                C.Ma=bsxfun(@times,C.Ma,Waux);
            end
            R.Ti.la=1./(abs(C.Ma).^2+0.01);%Inverse from 0.01 to 400-Direct from 0.0025 to 100
            A.Se=(sum(abs(conj(E.Sf).*E.Sf),4)+R.Ti.la).^(-1);       
            C=[];
        else
            if size(E.Sf,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction
                NMa=[E.Uf{1}.NX E.Uf{2}.NX];
                for n=1:2
                    if ~isempty(E.Uf{n});C.Ma=fold(C.Ma,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
                end           
                C.Ma=single(C.Ma>0);
                C.Ma=resampling(C.Ma,NMa,2);
            end
            if isfield(rec,'W');C.Ma=bsxfun(@times,C.Ma,single(E.WW>Alg.parE.eigSc));end
        end

        %Solver
        vDyn=[5:6 8 10:rec.Plan.NDims];%These are the dynamic directions
        vEch=[7 9];%These are the echo directions
        vExc=1:4;%Other directions
        y=permute(y,[vExc vEch vDyn]);
        if rec.Par.Mine.Modal==7;y=mean(y,12);end%Average dynamics
        NY=size(y);NY(end+1:rec.Plan.NDims)=1;
        y=reshape(y,[NY(1:4) prod(NY(5:6)) prod(NY(7:end))]);
        NX=size(E.Sf);NX(4:5)=[prod(NY(5:6)) prod(NY(7:end))];NX(6)=1;
        overDec=Alg.OverDec(1:3);
        overDec(overDec<0)=1;
        NX(1:3)=round(NX(1:3)./overDec);

        %MEMORY FOR RECONSTRUCTIONS, G-FACTORS AND CHI2-FACTORS
        x=single(zeros(NX));if gpuIn;x=gpuArray(x);end
        if Alg.EstimCFact
            c=single(zeros(NX));if gpuIn;c=gpuArray(c);end
        end
        if Alg.EstimGFact(1)
            if strcmp(Alg.solverType,'CG')
                if ~isfield(rec,'G');g=single(zeros(NX(1:3)));else g=[];end
            else
                g=single(zeros(NX));
            end
            if gpuIn;g=gpuArray(g);end
        end
        if isfield(rec,'G') && size(rec.G,4)==NPE && ~sepno;estGF(1)=0;end
        
        if NX(5)~=1 && Alg.parX.refineMask%RECONSTRUCTION OF AVERAGED INFORMATION---NOT WORKING FOR BLIPPED
            voxsiz=rec.Enc.AcqVoxelSize;
            voxsiz(3)=voxsiz(3)+rec.Par.Labels.SliceGaps(1);
            parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
            parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
            parS.Otsu=[];%Binary vector of components for multilevel extraction (it picks those components with 1)
            parS.nErode=0;%Erosion for masking (in mm)
            parS.nDilate=[64.1 2.1 1.1].*voxsiz;%Dilation for masking (in mm)
            parS.conComp=2;%Whether to get the largest connected component after erosion  
            
            if rec.Par.Mine.AdHocArray(8)~=11
                yen=multDimMea(y,5:6);E.cc=1;
                
                if strcmp(Alg.solverType,'CG');xaux=CGsolver(yen,E,EH,A,C,R,[],[],[],Alg.TolerSolve);
                elseif strcmp(Alg.solverType,'IRWLS');xaux=IRWLSsolver(yen,E,EH,Alg.IRWLSPower,4,1e-2,A,C,R,[],[],Alg.TolerSolve,[],[],1);
                end
                C.Ma=refineMask(xaux,parS,voxsiz);
                rec.M=C.Ma;
                rec.S=E.Sf;
            else                
                yen=y;
                if isfield(EH,'Ub')%Sense unfolding (third dimension)             
                    for n=3        
                        if isfield(EH,'Bb');yen=repmat(yen,[1 1 EH.Ub{n}.NX/EH.Ub{n}.NY]);elseif ~isempty(EH.Ub{n});yen=ifold(yen,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
                    end   
                end            
                %if isfield(EH,'Gh');yen=bsxfun(@times,yen,EH.Gh.Ab);end%Nyquist ghosting decoding
                if isfield(EH,'Bb');yen=bsxfun(@times,yen,EH.Bb);end%MB decoding
                yen=multDimMea(yen,5:6);                
                if isfield(EH,'Eb');yen=aplGPU(EH.Eb,yen,E.pe);end%Fourier 
                
                for rm=1:1            
                    ES=E;
                    ESH=EH;
                    ES=rmfield(ES,{'Bl','Ef','Bf'});
                    ESH=rmfield(ESH,{'Eb','Bb'});              
                    ES.Uf{3}=[];ES.UAf{3}=[];
                    ESH.Ub{3}=[];ESH.UAb{3}=[];
                    
                    if strcmp(Alg.solverType,'CG');xaux=CGsolver(yen,ES,ESH,A,C,R,[],[],[],Alg.TolerSolve);
                    else xaux=IRWLSsolver(yen,ES,ESH,Alg.IRWLSPower,4,1e-2,A,C,R,[],[],Alg.TolerSolve,[],[],1);
                    end
                    %xaux=bsxfun(@times,C.Ma,bsxfun(@times,A.Se,sum(bsxfun(@times,conj(E.Sf),yen),4)));
                    xaux=margosianFilter(xaux,rec.Enc,1);                  
                                   
                    recS=rec;               
                    recS.x=xaux;
                    recS.y=yen;
                    recS.Enc.RecSize=size(xaux);
                    recS.Enc.FOVSize=size(yen);             
                    recS=solveS(recS);                              
                
                    E.Sf=recS.S;recS=[];
                    A.Se=(sum(abs(conj(E.Sf).*E.Sf),4)+1e-9).^(-1);       
                    C.Ma=refineMask(xaux,parS,voxsiz);
                    rec.M=C.Ma;
                    rec.S=E.Sf;
                end
            end              
        end

        %dev=gpuDevice;
        for e=1:NX(4);E.ec=e;%Echo directions            
            for n=1:NX(5);if rec.Par.Mine.AdHocArray(8)==11;E.cc=n;else E.cc=1;end%Cardiac directions
                yen=dynInd(y,[e n],5:6);
                if isfield(rec.Par.Mine,'ShiftIsoc') && Alg.JointSplitScans
                    E.Gh.Pf=dynInd(P.PV,indSplit(n),12);EH.Gh.Pb=dynInd(P.PVH,indSplit(n),12);
                    E.Bf=dynInd(B.BV,indSplit(n),12);EH.Bb=dynInd(B.BVH,indSplit(n),12);
                end
                
                
                if all(multDimMax(yen,[1 2 4])~=0)
                    if Alg.PE==1;[yen,E,EH,A,C,R]=flipPE(yen,E,EH,A,C,R);end                     
                    if isfield(E,'Gh')
                        [E,EH]=computeEGhA(E,EH);
                        if size(EH.Gh.Ab,EH.pe)~=size(yen,EH.pe);fprintf('PE dimensions not matched between Nyquist ghost calibration and data for file %s\n',rec.Names.Name);rec.Fail=1;return;end
                        if size(EH.Gh.Ab,3-EH.pe)~=size(yen,3-EH.pe);fprintf('RO dimensions not matched between Nyquist ghost calibration and data for file %s\n',rec.Names.Name);rec.Fail=1;return;end                            
                    end
                    
                    if ~sepno && (e~=1 || n~=1);estGF(1)=0;end
                    %if (isfield(rec,'G') && ~isempty(rec.G)) || (e~=1 || n~=1);estGF(1)=0;end
                    
                    if Alg.artifRobust==1;[E,EH]=artifactEstimation(E,EH,yen,Alg.artifSigma);end
                    if strcmp(Alg.solverType,'CG');[xaux,~,gaux,caux,~,~]=CGsolver(yen,E,EH,A,C,R,[],[],[],Alg.TolerSolve,estGF*single(~unfat),Alg.EstimCFact);
                    elseif strcmp(Alg.solverType,'IRWLS');[xaux,~,gaux,caux,~]=IRWLSsolver(yen,E,EH,Alg.IRWLSPower,4,1e-2,A,C,R,[],[],Alg.TolerSolve,estGF*single(~unfat),Alg.EstimCFact);
                    else error('No valid solver (%s) has been selected',Alg.solverType);
                    end
                    
                    %if Alg.artifRobust==2
                    %    [E,EH]=artifactEstimation(E,EH,yen,Alg.artifSigma,caux);
                    %    if strcmp(Alg.solverType,'CG');[xaux,~,gaux,caux]=CGsolver(yen,E,EH,A,C,R,[],[],[],Alg.TolerSolve,estGF,Alg.EstimCFact);
                    %    elseif strcmp(Alg.solverType,'IRWLS');[xaux,~,gaux,caux]=IRWLSsolver(yen,E,EH,Alg.IRWLSPower,4,1e-2,A,C,R,[],[],Alg.TolerSolve,estGF,Alg.EstimCFact,1);
                    %    else error('No valid solver (%s) has been selected',Alg.solverType);
                    %    end
                    %end
                    
                    if unfat                                                                                       
                        laPrev=R.Ti.la;
                        SePrev=A.Se;
                        SfPrev=E.Sf;
                        
                        caux=sqrt(caux);
                        M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX],'tukey',[0.5 1],gpuIn,1,1);
                        caux=filtering(caux,M,1);                        
                        caux=caux-median(caux(:));
                        sigma=1.4826*median(abs(caux(:)));
                        tau=1.345*sigma;               

                        cauxa=caux/(5*tau)-1;%%%TEST AS IT IS AND COMPARE WITH PREVIOUS ONE IN TERMS OF ARTIFACT REMOVAL
                        cauxa=0.5+0.5*tanh(3*cauxa);%Standard sigmoid---1 indicates almost sure artifact, 0 indicates no artifact
                        cauxa=abs(filtering(cauxa,M,1));

                        Waux=(1-cauxa)*Alg.parE.eigSc(1)+cauxa*Alg.parE.eigSc(2);
                        Waux=2*bsxfun(@times,bsxfun(@minus,E.WW,Waux),1./(1-Waux))-1;
                        Waux=0.5+0.5*tanh(3*Waux);%"Standard" sigmoid                        

                        if Alg.UnfoldFat==2                      
                            zf=1;
                            if rec.Par.Labels.WFS==0;SH{E.pe}=-round(29.5507);else SH{E.pe}=-round(rec.Par.Labels.WFS);end%Water-Fat shift                                                      
                            E.Sf=cat(6,E.Sf,shifting(dynInd(E.Sf,1,6),SH,[],zf));                          
                            cauxb=caux/(20*tau)-1;%%%TEST AS IT IS AND COMPARE WITH PREVIOUS ONE IN TERMS OF ARTIFACT REMOVAL
                            cauxb=0.5+0.5*tanh(3*cauxb);%Standard sigmoid---1 indicates almost sure artifact, 0 indicates no artifact
                            cauxb=abs(filtering(cauxb,M,1));                      
                            Waux=cat(6,Waux,cauxb.*shifting(dynInd(Waux,1,6),SH,[],zf));             
                        end                                
                        Waux=Alg.parE.lambda*Waux;                        
                        R.Ti.la=1./(abs(Waux).^2+0.01);%Inverse from 0.01 to 400-Direct from 0.0025 to 100
                        A.Se=(sum(abs(conj(E.Sf).*E.Sf),4)+R.Ti.la).^(-1);                       
                        
                        if strcmp(Alg.solverType,'CG');[xaux,~,gaux,caux]=CGsolver(yen,E,EH,A,C,R,[],[],[],Alg.TolerSolve,estGF,Alg.EstimCFact);
                        elseif strcmp(Alg.solverType,'IRWLS');[xaux,~,gaux,caux]=IRWLSsolver(yen,E,EH,Alg.IRWLSPower,4,1e-2,A,C,R,[],[],Alg.TolerSolve,estGF,Alg.EstimCFact);
                        else error('No valid solver (%s) has been selected',Alg.solverType);   
                        end
                    end                    
                    if Alg.GhosCorRec && ismember(rec.Par.Mine.Modal,9:10);xaux=solveXP(xaux,yen,E,EH,A,C,R,[],[],Alg.TolerSolve);end
                    if unfat
                        A.Se=SePrev;
                        R.Ti.la=laPrev;
                        E.Sf=SfPrev;
                    end
                    
                    if ~isempty(xaux);xaux=dynInd(xaux,1,6);end%FOR ESPIRIT
                    if ~isempty(gaux);gaux=dynInd(gaux,1,6);end%FOR ESPIRIT

                    if Alg.PE==1;[~,E,EH,A,C,R,xaux,gaux,caux]=flipPE(yen,E,EH,A,C,R,xaux,gaux,caux);end
                    if Alg.parV.visReconstruction==1;visReconstruction(xaux);end                                       

                    if Alg.MargosianFilter%PARTIAL FOURIER   -BEFORE OVERDECODING DUE TO BOUNDARY ARTIFACTS IN VOLUMETRIC SCANS
                        Enc=rec.Enc;
                        xaux=margosianFilter(xaux,Enc,e); 
                    end                   
                    xaux=removeOverencoding(xaux,Alg.OverDec);%REMOVE OVERDECODING
                    if ~isempty(caux);caux=removeOverencoding(caux,Alg.OverDec);end%REMOVE OVERDECODING
                    x=dynInd(x,[e n],4:5,xaux);
                    if ~isempty(caux);c=dynInd(c,[e n],4:5,caux);end
                    if ~isempty(gaux)                        
                        if sepno || (e==1 && n==1 && (~isfield(rec,'G') || p>size(rec.G,4)));gaux=removeOverencoding(gaux,Alg.OverDec);end%REMOVE OVERDECODING, we assume PARTIAL FOURIER is not required
                        if sepno;g=dynInd(g,[e n],4:5,gaux);elseif e==1 && n==1 && (~isfield(rec,'G') || p>size(rec.G,4));g=gaux;end
                    end
                else
                    fprintf('Data looks corrupted for volume %d echo %d\n',n,e);
                end
            end
        end        
        E=[];EH=[];C=[];A=[];R=[];yen=[];
        x=resSub(x,[4 5]);
        x=dynInd(x,multDimMax(x,1:3)~=0,4);
        if Alg.EstimCFact
            c=resSub(c,[4 5]);
            c=dynInd(c,multDimMax(c,1:3)~=0,4);
        end
        if ~gpu;x=gather(x);end        
        if isempty(xp);xp=x;else xp=cat(4,xp,x);end%This differentiation is necessary in the gpu, probably because of a bug in cat
        if Alg.EstimCFact
            if isempty(cp);cp=c;else cp=cat(4,cp,c);end%This differentiation is necessary in the gpu, probably because of a bug in cat
        end
        if Alg.EstimGFact(1)
            if ~gpu;g=gather(g);end
            g=resSub(g,[4 5]);
            if isempty(gp);gp=g;else gp=cat(4,gp,g);end
        end
    end
end
if NPE>1
    xp=dynInd(xp,vertcat(indPEV{:}),4,xp);
    if Alg.EstimCFact;cp=dynInd(cp,vertcat(indPEV{:}),4,cp);end
end
if ~any(rec.Dyn.Typ2Rec==12);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,12);rec.Dyn.Typ2Wri(12)=1;end
if Alg.EstimCFact && ~any(rec.Dyn.Typ2Rec==10);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,10);rec.Dyn.Typ2Wri(10)=1;end
if rec.Par.Mine.Modal==9;xp=gather(xp);end%To save memory for long fMRI studies
if ~isfield(rec,'x') || isempty(rec.x);rec.x=xp;else rec.x=cat(4,rec.x,xp);end%This differentiation is necessary in the gpu, probably because of a bug in cat
if Alg.EstimCFact
    if rec.Par.Mine.Modal==9;cp=gather(cp);end%To save memory for long fMRI studies
    if ~isfield(rec,'X') || isempty(rec.X);rec.X=cp;else rec.X=cat(4,rec.X,cp);end%This differentiation is necessary in the gpu, probably because of a bug in cat
end
if Alg.EstimGFact(1) 
    if sepno
        if NPE>1;gp=dynInd(gp,vertcat(indPEV{:}),4,gp);end        
        if ~isfield(rec,'G') || isempty(rec.G);rec.G=gp;else rec.G=cat(4,rec.G,gp);end%This differentiation is necessary in the gpu, probably because of a bug in cat
    elseif ~isfield(rec,'G')
        rec.G=gp;
    elseif size(rec.G,4)~=NPE
        rec.G=cat(4,rec.G,gp);
    end
end
