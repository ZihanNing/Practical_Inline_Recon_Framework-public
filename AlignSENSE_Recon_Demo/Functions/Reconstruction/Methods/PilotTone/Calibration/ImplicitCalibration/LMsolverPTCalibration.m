
function [E,x,PT] = LMsolverPTCalibration(y,E,EH,x,C,deb,PT,parPT,nIt,tol)

%LMSOLVER Performs a Levenger-Marquardt iteration on a series of 
%   parameters of the encoding matrix to minimize a least squares 
%   backprojection of the reconstruction to the measured data
%   E=LMsolver(Y,E,X,{C},{DEB})
%   * Y is the measured data
%   * E is the encoding structure
%   * X is the reconstructed data
%   * {C} is a constrain structure
%   * {DEB} indicates whether to print information about convergence
%   ** E is the modified encoding structure
%

if nargin<5;C=[];end
if nargin<6 || isempty(deb);deb=2;end
if nargin<7 || isempty(PT);warning('LMsolverPTCalibration:: No PT structure detected.\nReturned original calibration parameters.');return;end
if nargin<8 || isempty(parPT);warning('LMsolverPTCalibration:: No parPT structure detected.\nReturned original calibration parameters.');return;end
if nargin<9 || isempty(nIt);nIt=10;end
if nargin<10 || isempty(tol);tol=1e-5;end

if isfield(EH,'We')&& ~isempty(EH.We) && parPT.Calibration.weightedCalibrationFlag
    We = EH.We; %EH.Mbe not included since already accomplished via E.nF region extraction
end
calibToDegrees=1;%Test for rot-ran imbalance

%% ========== TEMP STORAGE
% close all
% pTimeResGrOrig = PT.pTimeResGr;
% %pTimeResGrOrig = cat(1,real(pTimeResGrOrig),imag(pTimeResGrOrig));
% AbOrig = PT.Ab;
% nIt=15;
% tol=0;
% E = rmfield(E,'nF');
% E.tL = zerosL(E.tL);
%AbOrig = cat(2,AbOrig,AbOrig);
% %% ========== TEMP SIMULATION
% isOffset = size(AbOrig,2)>size(pTimeResGrOrig,1);
% rng(22);
% NChaTest=10;
% paramSelect=1:6; paramsDiscard=dynInd(1:6,~ismember(1:6,paramSelect),2);
% 
% PT.NChaPTRec = NChaTest;
% PT.pTimeResGr = dynInd(pTimeResGrOrig,1:NChaTest,1);
% 
% PT.Ab = dynInd(AbOrig,1:(NChaTest+isOffset),2);
% 
% PT.pTimeResGr = (PT.pTimeResGr).*exp(1i*0*pi/10*2*plugNoise(abs(PT.pTimeResGr),1,0));
% PT.pTimeResGr =  PT.pTimeResGr(1:NChaTest,:) + 1i* pTimeResGrOrig(33:33+NChaTest-1,:) ;
% isComplexFlag = ~isreal(PT.pTimeResGr);
% %PT.pTimeResGr = bsxfun(@minus,PT.pTimeResGr,multDimMea(PT.pTimeResGr,2));
% PT.pTimeResGr = PT.pTimeResGr + .01*plugNoise(PT.pTimeResGr,1,0);
% 
% PTGT=PT;
% amp = 1200;
% PTGT.Ab = amp*plugNoise((PT.Ab(:,1:NChaTest)),0,0)/NChaTest%.*exp(1i*0*pi/6*2*plugNoise(PT.Ab,1,0));
% PTGT.Ab = dynInd(PTGT.Ab, 4:6, 1, pi/180*dynInd(PTGT.Ab,4:6,1));
% %if isComplexFlag ;PTGT.Ab = PTGT.Ab .* ang
% PTGT.Ab = dynInd(PTGT.Ab, paramsDiscard, 1, 0);
% %if isOffset;PTGT.Ab = dynInd(PTGT.Ab, NChaTest+1, 2, dynInd(PTGT.Ab, NChaTest+1, 2)/5);end
% 
% %PTGT.Ab = ones(size(PT.Ab));
% TGT=pilotTonePrediction (PTGT.Ab, PT.pTimeResGr, 'PT2T','backward', PT.NChaPTRec, size(E.Tr,6));
% E.Tr=TGT; 
% EGT=E; 
% y = encode(x,EGT);
% if isequal(paramSelect,1:3)
%     PT.Ab = PTGT.Ab + amp/10*[ ones([1 3]) 0*ones([1 3])]'.*rand([NChaTest+isOffset 6])';
% elseif isequal(paramSelect,4:6)
%     PT.Ab = PTGT.Ab + amp/5*[ 0*ones([1 3]) pi/180*ones([1 3])]'.*rand([NChaTest+isOffset 6])';
% else
%     PT.Ab = PTGT.Ab + amp/50*[ ones([1 3]) pi/180*ones([1 3])]'.*rand([NChaTest+isOffset 6])';
% end
% if isOffset;PT.Ab = dynInd(PT.Ab, NChaTest+1, 2, dynInd(PT.Ab, NChaTest+1, 2)/5);end
% 
% E.Tr = pilotTonePrediction (PT.Ab, PT.pTimeResGr, 'PT2T','backward', PT.NChaPTRec, size(E.Tr,6));
% visMotion(TGT,[],[],0,[],[],[],[],[],65);
% visMotion(TGT-E.Tr,[],[],0,[],[],[],[],[],66);
% %linePlotA(y, E, EH, x, PT.Ab, PTGT.Ab, [], [], 100, 11, PT.pTimeResGr );

%% ==================

if isfield(E,'Tr')              
    %GENERAL PARAMETERS AND ARRAYS
    gpu=isa(x,'gpuArray');
    T=E.Tr;E=rmfield(E,'Tr'); %YB: need to remove it otherwise the encode will perform the tranfsormation for a second time
    if isfield(E,'Dl');TD=E.Dl.Tr;end
    if isfield(E,'Ds');TD=E.Ds.Tr;end
    
    if isfield(E,'Db'); Db = E.Db; E = rmfield(E,'Db');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dc'); Dc = E.Dc; E = rmfield(E,'Dc');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dl'); Dl = E.Dl; E = rmfield(E,'Dl');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Ds'); Ds = E.Ds; E = rmfield(E,'Ds');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'B0Exp'); B0Exp = E.B0Exp; E = rmfield(E,'B0Exp');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dbs'); Dbs = E.Dbs; E = rmfield(E,'Dbs');end%Could include this part in encode.m, but it will be skipped there as it in in the if E.Tr....end
    if isfield(E,'B1m'); B1m = E.B1m; E = rmfield(E,'B1m');end
    if isfield(E,'Shim'); Shim = E.Shim; E = rmfield(E,'Shim');end  
    
    NT=size(T);
    NX=size(x);    
    dimG=max(6,length(NT));dimS=find(NT(1:dimG-1)~=1);if isempty(dimS);dimS=dimG-1;end%Dimensions of motion states and parameters of the transform 
    ndS=NT(dimS);ndG=NT(dimG);NTS=length(dimS);
    ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
    ndH=size(ha,2);   
    %flagw=zeros(NT(1:dimS(end)));%PTHandling: This was for LMsolver
    flagw=0;%PTHandling - only a single energy since not decoupled anymore
    En=single(zeros([ndS 1]));EnPrev=En;EnPrevF=EnPrev;EnF=En;   
    dH=single(zeros([ndS ndH]));dG=single(zeros([ndS ndG]));dGEff=dG; %YB: For every state you have 21 (6+5+4+3+2+1) Hessian elements and 6 gradient terms  
    multA=1.2;multB=2;%YB: Factors to divide/multiply the weight that regularizes the Hessian matrix when respectively E(i)<E(i-1) and  E(i)>E(i-1)    
    NElY=numel(y);
    if isfield(E,'nF');E.Sf=dynInd(E.Sf,E.nF,3);y=dynInd(y,E.nF,3);end%Extract ROI in the readout direction
    if isfield(E,'Ps') && ~isempty(E.Ps);E.Sf=bsxfun(@times,E.Sf,E.Ps);y=bsxfun(@times,y,E.Ps);end%Preconditioner of the coils, deprecated
    
    %PTHandling - convert PT signal and calbiration matrix in new matrix notation so we can use linear algebra equations
    A=[];indA=[];if exist('PTGT','var');AGT=[];end
    isOffset = size(PT.Ab,2)>PT.NChaPTRec; %Offset in calibration;
    p=zeros([ndG ndG*(PT.NChaPTRec+isOffset) ndS], 'like', PT.pTimeResGr);
        
    for i=1:ndG
        A=cat(1,A,PT.Ab(i,:).');
        if exist('AGT','var');AGT=cat(1,AGT,PTGT.Ab(i,:).');end
        indA = cat(1, indA, i*ones([length(PT.Ab(i,:)) 1]));
    end
    for i=1:ndG
        for ii=1:ndS
            if isOffset %Offset in calibration
                p = dynInd(p,{i,indA==i, ii},1:3,cat(2,dynInd(PT.pTimeResGrComp,{ii},2).',1));
            else
                p = dynInd(p,{i,indA==i, ii},1:3,dynInd(PT.pTimeResGrComp,{ii},2).');
            end
        end %Inefficient!! Change to more efficient implementation
    end
    %[A,p,indA]=reshapeP(PT.Ab,PT.pTimeResGrComp,'mat2vec'); 
    %if exist('PTGT','var');AGT=reshapeP(PTGT.Ab,[],'mat2vec');end
    E.wCalib = E.winit*onesL(flagw);%PTHandling: for calibration - a scalar
    E.w = E.winit*onesL(E.w);%reset - was not in LMsolver
    E.cT = zerosL(E.cT);%PTHandling: after multiple LMsolver calls, there might be a lot converged, so reset
    disp('convergence reset ---------- revieeeeeeeeeeeeeeeeeeeew')

    %BLOCK SIZES AND CONVERGENCE VALUES FOR THE TRANSFORM
    ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;ET.kG=E.kG;ET.rkG=E.rkG;   
    if isfield(E,'Bm');ET.Bm=E.Bm;ET.bS(1)=1;end%YB: in encode.m = slice mask
    if isfield(E,'Fms');ET.Sf=E.Sf;ET.Fms=E.Fms;ET.bS(1)=NX(3);E.dS(1)=1;end%YB: in encode.m = Fourier encoding multislice
    if isfield(E,'ZSl');ET.ZSl=E.ZSl;E.ZSl=-E.ZSl;end%YB: in encode.m = Slab extraction
    ET.cTFull=E.cT;       

    for n=1:nIt
        %nResetImplicit = 5;%PTHandling: to reset the number of shots to use for gradient calculation
        %if mod(n,nResetImplicit)==0;ET.cTFull = zerosL(ET.cTFull);end%PTHandling: in previous iteration,ET.cTFull might be set to 1 if some states are not updated a lot
        ET.cTFull = zerosL(ET.cTFull);
        fprintf('LMsolverPTCalibration: Iteration %d\n',n)
        fprintf('States used for gradient calculations: %d/%d\n',multDimSum(~ET.cTFull),length(ET.cTFull))
        
        flagw=0;%PTHandling: In LMsolver there was only 1 iteration, so it was not explicitely repeated
        if NTS==1;NB=1;else NB=NT(dimS(1));end
        for b=1:NB
            if NB==1;ET.cT=ET.cTFull;else ET.cT=dynInd(ET.cTFull,b,dimS(1));end
            if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

            %COMPUTE THE JACOBIAN AND STORE VALUES
            ind2EstOr=find(~ET.cT(:));NEst=length(ind2EstOr); %YB: indices which state to update (so not converged)
            if isfield(E,'Fs');ry=zeros([size(y,1) 1],'like',real(y));else; ry=zeros([NT(6) 1],'like',real(y));end

            if exist('Ds','var') && isfield(Ds,'B0')% Susceptibility model with voxel basis
                x = bsxfun(@times, x , exp(1i*2*pi*Ds.TE*Ds.B0));
            end

            for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2EstOr(vA);
                if isfield(E,'Fms')
                   [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
                   E.Sf=dynInd(ET.Sf,vA,3);            
                end                               
                if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
                if NB==1;vT=vA;else vT={b,vA};end
                Ti=dynInd(T,vT,dimS);
                [Tf,Tg]=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,1,1,1);%Transform factors 

                xT=x;
                if exist('B0Exp','var');xT=bsxfun(@times,xT, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end  
                if exist('Db','var');xT=bsxfun(@times,xT,dephaseBasis( Db.B, dynInd(Db.cr,vA(vA<=E.NMs),5), size(xT),Db.TE));end
                if exist('Dc','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
                %if exist('Dl','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dl.d,6),Dl.D,Dl.TE));end   % Linear model with voxel model
                if exist('Dl','var')
                    TiD = dynInd(TD,vT,dimS);
                    xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD,Dl.f,Dl.D,Dl.TE));
                end
                if exist('Ds','var') && isfield(Ds,'D')
                    TiD=dynInd(TD,vT,dimS);
                    xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Ds.f,Ds.D,Ds.TE));
                end   % Susceptibility model with voxel basis
                if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xT=bsxfun(@times,xT,intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xT),Db,Dl));end

                if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
                if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                       
                [xT,xTG]=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);

                if exist('Dbs','var');xT=bsxfun(@times,xT,dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xT),Dbs.TE));end
                if exist('B1m','var');xT=bsxfun(@times,xT,exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xT),[],1)));end   
                if exist('Shim','var');xT=bsxfun(@times,xT,exp(+1i*2*pi*Shim.TE *Shim.B0) );end   

                if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end%YB:Taking out E.nF only after transformation, otherwise transformations not consistent with rest of algorithm
                if isfield(E,'Fs')
                    E.vA=vA;
                    E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                    xT=encode(xT,E)-dynInd(y,indY,1); %YB:This is to create the term w in equation 11 of Lucilio's first paper
                elseif isfield(E,'Bm')
                    E.Bm=dynInd(ET.Bm,vA,6);
                    xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
                else
                    xT=encode(xT,E)-dynInd(y,vT,[5 3]);
                end
                
                if exist('We','var'); xT = bsxfun(@times,xT,dynInd(We,indY,1) );end%PTHandling: weighting parts of k-space

                if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
                if isfield(E,'Fs')
                    ry(indY)=normm(xT,[],3:4)/NElY;
                    EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
                else
                    ry(vA)=normm(xT)/NElY;
                    EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(ry(vA)));
                end
                if isfield(E,'Fd') && ~isempty(E.Fd)%Filtering
                    xT=fftGPU(xT,3);
                    xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                    EnPrev=dynInd(EnPrev,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
                else
                    EnPrev=dynInd(EnPrev,vT,1:NTS,dynInd(EnPrevF,vT,1:NTS));
                end

                xTG=sincRigidTransformGradient(xTG,Tf,Tg,E.Fof,E.Fob); % This creates outer 2 operations of w_l in Lucilios 1st paper (still needs encoding!)

                for g=1:ndG %YB: Building Gradient (with the use of the Jacobian)
                    %Apply phase in scanner reference to Gradient as well
                    if exist('Dbs','var');xTG{g}=bsxfun(@times,xTG{g},dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xTG{g}),Dbs.TE));end
                    if exist('B1m','var');xTG{g}=bsxfun(@times,xTG{g},exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xTG{g}),[],1)));end   
                    if exist('Shim','var');xTG{g}=bsxfun(@times,xTG{g},exp(+1i*2*pi*Shim.TE *Shim.B0) );end 

                    if isfield(E,'nF');xTG{g}=dynInd(xTG{g},E.nF,3);end
                    xTG{g}=encode(xTG{g},E); % YB: This is the encoding step still needed to become term w_l in 11th equation in 1st paper
                    if isfield(E,'Pf') && ~isempty(E.Pf);xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Pf,indY,1));end  %Precondition using all the residuals, deprecated
                    if isfield(E,'Fd') && ~isempty(E.Fd)
                        xTG{g}=fftGPU(xTG{g},3);
                        xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Fd,indY,1));
                    end

                    if isfield(E,'Fs')
                        dG(vA,g)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{g}.*conj(xT)),3:4)),vA,1));%YB: This is the gradient calculation in Eq. 11 in the alignedSENSE paper
                    elseif isfield(E,'Bm')
                        dG(vA,g)=gather(multDimSum(real(xTG{g}.*conj(xT)),1:6));
                    else
                        dG=dynInd(dG,[vT g],1:NTS+1,gather(permute(multDimSum(real(xTG{g}.*conj(xT)),[1:2 4]),[5 3 1 2 4])));
                    end
                end
                for h=1:ndH %YB: Building Hessian - approximated by JHJ = different from AlignedRecon paper
                    if isfield(E,'Fs')
                        dH(vA,h)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),3:4)),vA,1));
                    elseif isfield(E,'Bm')
                        dH(vA,h)=gather(multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),1:6));
                    else
                        dH=dynInd(dH,[vT h],1:NTS+1,gather(permute(multDimSum(real(xTG{ha(1,h)}.*conj(xTG{ha(2,h)})),[1:2 4]),[5 3 1 2 4])));
                    end
                end
            end;xTG=[];
        end

        %UPDATE
        En=EnPrev;EnF=EnPrevF;%Initialise and fill with values from new update
        MH=single(eye(ndG));% YB: This is the actual matrix of the Hessian with symmetry + the LM parameters included on the diagonal
        EnPrevOrig = EnPrev;%PTHandling: need to store previous energy
        fina=0;
        while ~fina   
            %BUILD HESSIAN MATRIX AND POTENTIAL UPDATE
            ry(:)=0;
            MHAll = cell([1 ndS]);%cell([1 multDimSum(~ET.cTFull(:) & E.w(:)<1e10)' ]);%PTHandling
            GJAll = cell([1 ndS]);%cell([1 multDimSum(~ET.cTFull(:) & E.w(:)<1e10)' ]);%PTHandling
            for a=find(~ET.cTFull(:) & E.w(:)<1e10)'            
                aa=ind2subV(NT(dimS),a);
                for h=1:ndH
                    if ha(1,h)==ha(2,h)
                        MH(ha(1,h),ha(2,h))=(1+E.wCalib)*dynInd(dH,[aa h],1:NTS+1);%PTHandling: only a single weight E.w for the calibration estimation
                        %MH(ha(1,h),ha(2,h))=(1+E.w(a))*dynInd(dH,[aa h],1:NTS+1);% YB: Here, the LM parameter w is already included, so no need for identity matrix    
                    else
                        MH(ha(1,h),ha(2,h))=dynInd(dH,[aa h],1:NTS+1);
                        MH(ha(2,h),ha(1,h))=dynInd(dH,[aa h],1:NTS+1);%YB:Hessian is symmetric
                    end
                end  
                MHAll{a} =  MH;%PTHandling
                GJAll{a} =  double(dynInd(dG,aa,1:NTS));%PTHandling
                if calibToDegrees
                    for dd=1:2;MHAll{a} = dynInd(MHAll{a},4:6,dd, convertRotation( dynInd(MHAll{a},4:6,dd),'deg','rad') );end
                    GJAll{a} = dynInd(GJAll{a},4:6,2,convertRotation( dynInd(GJAll{a},4:6,2),'deg','rad'));
                end
                %GJ=single(double(MH)\GJ(:));% YB: This is where you take the inverse of the Hessian
                %GJ=resPop(GJ,1,[],NTS+1);
                %dGEff=dynInd(dGEff,aa,1:NTS,-(E.winit/E.w(a))*GJ);%YB: this is initialised with 0, so converged shots will have Tupr=0
            end

            %PTHandling: sum over gradient and hesssian - divide by number of segments to allow for only computing both on a subset of NStates
            for i=find(~cellfun(@isempty,GJAll)); if ~isempty(GJAll{i});GJAll{i}=  dynInd(p,i,3)'*GJAll{i}.';end;end 
            GJCalib = (multDimMea(cat(4, GJAll{:}),4));GJAll=[];

            for i=find(~cellfun(@isempty,MHAll)); if~isempty(MHAll{i}); MHAll{i}= dynInd(p,i,3)'*MHAll{i}'*dynInd(p,i,3);end; end 
            MHCalib = (multDimMea(cat(4, MHAll{:}),4));MHAll=[];

            Aupr= -(E.winit/E.wCalib)*single(double(MHCalib)\GJCalib(:));
            %Aupr = real(Aupr);
            if isempty(Aupr)
               fprintf('---------------------------------------------------------------------------------------Aupr empty.\n');
               Aup=A;AupNew = zerosL(PT.Ab);for i=1:ndG; AupNew = dynInd(AupNew,i,1, dynInd(gather(Aup),indA==i,1).');end 
               Tup = pilotTonePrediction (AupNew, PT.pTimeResGrComp, 'PT2T','backward', PT.NChaPTRec, ndG);
               if exist('Dl','var');[~, TupD] = transformationConversion([], Tup, Dl);end
               if exist('Ds','var') && isfield(Ds,'D');[~, TupD] = transformationConversion([], Tup, Ds);end
               if ~exist('EnList','var');EnList=[];end
               break
            end
            Aup = A + Aupr;
            %AupNew = reshapeP(Aup,[],'vec2mat',indA);
            AupNew = zerosL(PT.Ab);for i=1:ndG; AupNew = dynInd(AupNew,i,1, dynInd(gather(Aup),indA==i,1).');end 
            AuprNew = zerosL(PT.Ab);for i=1:ndG; AuprNew = dynInd(AuprNew,i,1, dynInd(gather(Aupr),indA==i,1).');end 
            Tup = pilotTonePrediction (AupNew, PT.pTimeResGrComp, 'PT2T','backward', PT.NChaPTRec, ndG);
            Tupr = pilotTonePrediction (AuprNew, PT.pTimeResGrComp, 'PT2T','backward', PT.NChaPTRec, ndG);
            %Tupr=shiftdim(dGEff,-(dimS(1)-1));
            %Tup=restrictTransform(Tup);%Rotation parameters between -pi and pi
            if exist('Dl','var');[~, TupD] = transformationConversion([], Tup, Dl);end
            if exist('Ds','var') && isfield(Ds,'D');[~, TupD] = transformationConversion([], Tup, Ds);end

            %CHECK ENERGY REDUCTION
            E.cTFull=(ET.cTFull | flagw);%If flagw is a scalar, it will never be 1 here (only useful when flagw is a vector)
            E.cTFull = zerosL(E.cTFull);%PTHandling: Need to compute energy on all shots (otherwise you might reduce 5 shots of interest and worsen 100 others)
            for b=1:NB
                if NB==1;E.cT=E.cTFull;else E.cT=dynInd(E.cTFull,b,dimS(1));end
                if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

                ind2Est=find(~E.cT(:));NEst=length(ind2Est);
                for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2Est(vA);
                    if isfield(E,'Fms')
                        [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
                        E.Sf=dynInd(ET.Sf,vA,3);                   
                    end                
                    if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
                    if NB==1;vT=vA;else vT={b,vA};end
                    Ti=dynInd(Tup,vT,dimS);
                    Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,0,1,1);

                    xT=x;  %For Ds.B0, already computed early on     
                    if exist('B0Exp','var');xT=bsxfun(@times,xT, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end
                    if exist('Db','var');xT=bsxfun(@times,xT,dephaseBasis( Db.B, dynInd(Db.cr,vA(vA<=E.NMs),5), NX,Db.TE));end
                    if exist('Dc','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
                    %if exist('Dl','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dl.d,6),Dl.D,Dl.TE));end % Linear model with voxel model
                    if exist('Dl','var')
                        TiD = dynInd(TupD,vT,dimS);
                        xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD,Dl.f,Dl.D,Dl.TE));
                    end
                    if exist('Ds','var') && isfield(Ds,'D') 
                        TiD = dynInd(TD,vT,dimS);
                        xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD,Ds.f,Ds.D,Ds.TE));
                    end   % Susceptibility model with voxel basis                
                    %if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xT=bsxfun(@times,xT,intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xT),Db,Dl));end

                    if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
                    if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                
                    xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);

                    if exist('Dbs','var');xT=bsxfun(@times,xT,dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xT),Dbs.TE));end
                    if exist('B1m','var');xT=bsxfun(@times,xT,exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xT),[],1)));end   
                    if exist('Shim','var');xT=bsxfun(@times,xT,exp(+1i*2*pi*Shim.TE *Shim.B0) );end   

                    if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end            
                    if isfield(E,'Fs')
                        E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                        E.vA=vA;
                        xT=encode(xT,E)-dynInd(y,indY,1);
                    elseif isfield(E,'Bm')
                        E.Bm=dynInd(ET.Bm,vA,6);
                        xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
                    else
                        xT=encode(xT,E)-dynInd(y,vT,[5 3]);
                    end
                    
                    if exist('We','var'); xT = bsxfun(@times,xT,dynInd(We,indY,1) );end

                    if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
                    if isfield(E,'Fs')
                        ry(indY)=normm(xT,[],3:4)/NElY;
                        EnF=dynInd(EnF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
                    else
                        ry(vA)=normm(xT)/NElY;
                        EnF=dynInd(EnF,vT,1:NTS,gather(ry(vA)));                                 
                    end
                    if isfield(E,'Fd') && ~isempty(E.Fd)
                        xT=fftGPU(xT,3);
                        xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                        En=dynInd(En,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
                    else
                        En=dynInd(En,vT,1:NTS,dynInd(EnF,vT,1:NTS));
                    end   
                end
            end
            
            EnAll = multDimMea(En,1);%PTHandling - not separable anymore - mean in case different number of shots considered
            EnPrevAll = multDimMea(EnPrevOrig,1);%PTHandling - not separable anymore
            En = EnAll*onesL(En);EnPrev = EnPrevAll*onesL(EnPrev);
            %En(E.w(:)>1e10)=EnPrev(E.w(:)>1e10);
            flagw = (EnAll<=EnPrevAll); % YB: flags whether energy reduction in that state, when this is all ones, you quit (fina=1)
           
            if deb>=2;fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(EnPrevF(:)),sum(EnF(:)));end  
            if any(~flagw(:)) %YB: not all converged yet
               E.w(En(:)>EnPrev(:) & ~ET.cTFull(:))=E.w(En(:)>EnPrev(:) & ~ET.cTFull(:))*multB;
               E.wCalib = E.wCalib*multB;
            else
                E.w(~ET.cTFull(:))=E.w(~ET.cTFull(:))/multA;
                E.w(E.w<1e-8)=multA*E.w(E.w<1e-8);%To avoid numeric instabilities 
                E.wCalib = E.wCalib/multA;
                E.wCalib(E.wCalib<1e-8)=multA*E.wCalib(E.wCalib<1e-8);%To avoid numeric instabilities 
                A=Aup;
                T=Tup;
                if exist('Dl','var');TD = TupD; end%YB
                if exist('Ds','var') && isfield(Ds,'D');TD = TupD; end
                fina=1;
                traDiff=abs(dynInd(Tupr,1:3,dimG));
                rotDiff=abs(convertRotation(dynInd(Tupr,4:6,dimG),'rad','deg'));  
                traDiffMax=multDimMax(traDiff,dimS);rotDiffMax=multDimMax(rotDiff,dimS);    
                ET.cTFull=zerosL(ET.cTFull);%PTHandling: disable ET.cTFull to accumulates across iterations, so possibly 
                ET.cTFull(max(traDiff,[],dimG)<E.tL(1) & max(rotDiff,[],dimG)<E.tL(2))=1; %YB: at this point you update the convergence of the motion states
                if deb>=2
                    fprintf('Maximum change in translation (vox): %s/ rotation (deg): %s\n',sprintf('%0.3f ',traDiffMax(:)),sprintf('%0.3f ',rotDiffMax(:)));
                    fprintf('Not converged motion states: %s of %s\n',num2str(prod(ndS)-sum(single(ET.cTFull(:)))),num2str(prod(ndS)));
                end              
            end
        end
        
        if deb>=2
            if ~exist('EnList','var');EnList=[];end
            EnList=cat(2,EnList,EnAll);
            h=figure(89); 
            set(h,'color','w');
            plot(1:length(EnList), EnList);
            xlabel('Number of iterations [#]');
            ylabel('Residuals [a.u.]');
            title('Convergenve plot for LMsolverPTCalibration');
        end
        if exist('AGT','var');visMotion(TGT-T,[],[],0,[],[],[],[],[],67);end
% AGTtt=AGT;AGTtt = zerosL(PT.Ab);for i=1:ndG; AGTtt = dynInd(AGTtt,i,1, dynInd(gather(AGT),indA==i,1).');end 
% size(PT.Ab)
% 
% figure;imshow(real(AGTtt),[])
% figure;imshow(real(PT.Ab),[])

        %PTHandling: Check convergence to abort iterations of LM update
        if all(ET.cTFull) || isempty(Aupr); break;end
        if exist('AGT','var');fprintf('%d\n',multDimMin(    abs((A-AGT)./(AGT+eps('single')) )  )); end
        E.cT(:)=0; disp('convergence reset ---------- revieeeeeeeeeeeeeeeeeeeew')
    end

    %PTHandling: convert PT signal and calbiration matrix back to original notation
    for i=1:ndG; PT.Ab = dynInd(PT.Ab,i,1, dynInd(gather(Aup),indA==i,1).');end 
        
    T=constrain(T,C);
    E.Tr=T;E.cT=ET.cTFull;E.bS=ET.bS;E.dS=ET.dS;E.oS=ET.oS;E.kG=ET.kG;E.rkG=ET.rkG;
    E.cT(:)=0; disp('convergence reset ---------- revieeeeeeeeeeeeeeeeeeeew')
    if exist('Dl','var');Dl.Tr = TD;end%Updated one
    if exist('Ds','var') && isfield(Ds, 'D');Ds.Tr = TD;end%Updated one
    if exist('B0Exp','var'); E.B0Exp = B0Exp;end%YB add again
    if exist('Db','var'); E.Db = Db;end%YB add again
    if exist('Dc','var'); E.Dc = Dc;end%YB add again
    if exist('Dl','var'); E.Dl = Dl;end%YB add again
    if exist('Ds','var'); E.Ds = Ds;end%YB add again
    if exist('Dbs','var'); E.Dbs = Dbs;end%YB add again
    if exist('B1m','var'); E.B1m = B1m;end%YB add again
    if exist('Shim','var'); E.Shim = Shim;end%YB add again

    if isfield(E,'Fms');E.Sf=ET.Sf;E.Fms=ET.Fms;end
    if isfield(E,'ZSl');E.ZSl=ET.ZSl;end
    if isfield(E,'Bm');E.Bm=ET.Bm;end
    if isfield(E,'vA');E=rmfield(E,'vA');end
    if isfield(E,'Fd');E=rmfield(E,'Fd');end   
    Enmin=min(EnPrev,En);
    if ~isfield(E,'EnCalib');E.EnCalib=Enmin;else E.EnCalib(ind2EstOr)=Enmin(ind2EstOr);end%To keep a record of the achieved energy
    %EH.Mbe = dynInd(EH.Mbe,E.nF,3);
    %linePlotA(y, E, EH, dynInd(x,E.nF,3), PT.Ab, PTGT.Ab, [], [], 100, 11, PT.pTimeResGr );
end
end
