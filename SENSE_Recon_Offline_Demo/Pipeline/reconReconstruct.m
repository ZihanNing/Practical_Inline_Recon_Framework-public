function rec=reconReconstruct(rec,typ,field)

%RECONRECONSTRUCT   Reconstructs data
%   REC=RECONINVERT(REC,TYP,{FIELD})
%   * REC is a recon structure
%   * TYP is the type of data to reconstruct, one of the following, 1->body, 
%2->surface, 3->data
%   * {FIELD} is the field with header info, it defaults to 'image', other
%   options are 'refscan'
%   ** REC is a recon structure
%

if nargin<3 || isempty(field);field='image';end

ND=16;
typV=['B' 'S' 'x'];
t=typV(typ);

if typ==1
    %STANDARDIZE DATA
    rec.(t)=standardizeCoils(rec.(t),rec.N.(t));
    %COMPRESS BODY
    rec.B=compressCoils(rec.B,1);
elseif typ==2
    %MAPPING FROM BODY COIL TO SURFACE COIL FOV
    
    %ESTIMATE SENSITIVITIES AND SOFT MASK
    %if size(rec.S,11)==2
    %    rec.PS=sqrt(mean(dynInd(rec.S,2,11).*conj(dynInd(rec.S,1,11)),4));
    %    rec.S=rec.S.*cat(11,signz(rec.PS),signz(conj(rec.PS)));        
    %end
    rec.xS=rec.S;
    if rec.Alg.averagePolaritySensitivities==2;rec.S=mean(rec.S,11);end
    if ~isfield(rec,'B');estB=1;else estB=0;end
    rec.B=[];rec.M=[];rec.W=[];
    size(rec.S)
    rec.S=dynInd(rec.S,1,8);%We extract first echo
    for s=1:size(rec.S,11)
        S=dynInd(rec.S,s,11);
        if estB
            if strcmp(rec.Alg.bodyCoilEstimationMethod,'RSoS');B=sqrt(normm(S,[],4));
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'CC');B=compressCoils(S,1);
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'BCC');B=blockCompressCoils(S,1);
            elseif strcmp(rec.Alg.bodyCoilEstimationMethod,'CCPh');B=sqrt(normm(S,[],4)).*signz(compressCoils(S,1));
            else error('Virtual body coil estimation method %s not contemplated',rec.Alg.bodyCoilEstimationMethod);
            end
        end
        if strcmp(rec.Alg.sensitivityEstimationMethod,'Filter')
            S=gather(S);B=gather(B);AcqDelta=rec.Enc.S.AcqDelta;maskThS=rec.Alg.maskThS;
            %save('/mnt/nfs/home/lcordero/Matlab/Data/dat.mat','S','B','AcqDelta','maskThS');
            %exit
            [S,M]=solveSFilt(B,S,rec.Enc.S.AcqDelta,rec.Alg.maskThS);
            
        elseif strcmp(rec.Alg.sensitivityEstimationMethod,'Standard')
            rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta;rph.Alg.parS=rec.Alg.parS;
            rph=solveS(rph);
            S=rph.S;M=rph.M;rph=[];
        elseif strcmp(rec.Alg.sensitivityEstimationMethod,'ESPIRIT')
            rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta;rph.Alg.parE=rec.Alg.parE;
            %rph.y=rph.y.*signz(B).*conj(signz(dynInd(rph.y,1,4)));%Phase correction
            %visSegment(angle(dynInd(rph.y,1,4)))
            rph=solveESPIRIT(rph);
            S=rph.S;M=single(rph.W>rec.Alg.parE.eigSc(1));W=rph.Wfull;rph=[];
            rec.W=cat(11,rec.W,W);
        else
            error('Sensitivity estimation method %s not contemplated',rec.Alg.sensitivityEstimationMethod);
        end
        rec.S=dynInd(rec.S,s,11,S);
        rec.B=cat(11,rec.B,B);
        rec.M=cat(11,rec.M,M);
    end
    if rec.Alg.averagePolaritySensitivities==1;rec.S=mean(rec.S,11);rec.B=mean(rec.B,11);end
    rec.M=max(rec.M,[],11);
    if useGPU;rec.S=gpuArray(rec.S);rec.xS=gpuArray(rec.xS);end
    rec.xS=sum(conj(rec.S).*rec.xS,4)./normm(rec.S,[],4);  
    rec.xS=resPop(rec.xS,11,[],5);rec.S=resPop(rec.S,11,[],5);rec.B=resPop(rec.B,11,[],5);
    if ~isempty(rec.W);rec.W=resPop(rec.W,11,[],5);end
    if size(rec.xS,5)>1
        rec.PS=sqrt(dynInd(rec.xS,2,5).*conj(dynInd(rec.xS,1,5)));rec.PS=gather(rec.PS);
        rec.Geom.PS.MS=rec.Geom.S.MS;rec.Geom.PS.MT=rec.Geom.S.MT;
    end
    rec.S=gather(rec.S);rec.xS=gather(rec.xS);
elseif typ==3 && strcmp(field,'refscan')
    %if useGPU;rec.S=gpuArray(rec.S);rec.x=gpuArray(rec.x);rec.N.x=gpuArray(rec.N.x);end
    %if ~isempty(rec.N.x)
    %    [rec.x,n]=standardizeCoils(rec.x,rec.N.x);
    %    rec.S=standardizeCoils(rec.S,n);
    %end
    if useGPU;rec.S=gpuArray(rec.S);rec.x=gpuArray(rec.x);end
    rec.Sx=mapVolume(rec.S,rec.x,rec.Geom.S.MT,rec.Geom.x.MT,[],[],0,'linear');rec.S=gather(rec.S);
    rec.xS=sum(conj(rec.Sx).*rec.x,4)./normm(rec.Sx,[],4);
    rec.xS=resPop(rec.xS,11,[],5);
    rec.PS=sqrt(dynInd(rec.xS,2,5).*conj(dynInd(rec.xS,1,5)));
    rec.xS=gather(rec.xS);rec.PS=gather(rec.PS);
    rec.Geom.PS.MS=rec.Geom.x.MS;rec.Geom.PS.MT=rec.Geom.x.MT;
    rec=rmfield(rec,'x');
else
    if useGPU;rec.S=gpuArray(rec.S);rec.Ay=gpuArray(rec.Ay);rec.Az=gpuArray(rec.Az);rec.x=gpuArray(rec.x);rec.N.x=gpuArray(rec.N.x);end
    if ~isempty(rec.N.x)
        [rec.x,n]=standardizeCoils(rec.x,rec.N.x);
        rec.S=standardizeCoils(rec.S,n);
    end
    NP=gather(max(rec.Ay(:)));%Number of polarities

    rec.Sx=mapVolume(rec.S,rec.x,rec.Geom.S.MT,rec.Geom.x.MT,[],[],0,'linear');rec.S=gather(rec.S);    
    if isfield(rec,'PS')
        NPS=size(rec.PS,1:2);
        H=buildFilter(NPS,'tukey',0.125,useGPU,1);
        rec.PS=filtering(rec.PS,H);
        NPS=size(rec.PS,1:3);
        H=buildFilter(2*NPS,'tukey',0.125,useGPU,1,1);
        rec.PS=filtering(rec.PS,H,1);
        rec.PSx=mapVolume(rec.PS,rec.x,rec.Geom.PS.MT,rec.Geom.x.MT,[],[],0,'linear');rec.PS=gather(rec.PS);
        rec.PSx(isnan(rec.PSx(:)))=1;
        rec.PSx=signz(rec.PSx);
        rec.Sx=rec.Sx.*flip(cat(5,rec.PSx,conj(rec.PSx)),5);
    else
        rec.Sx=repmat(rec.Sx,[1 1 1 1 NP]);
    end
    if isfield(rec,'xS')
        if useGPU;rec.xS=gpuArray(rec.xS);end
        parS=rec.Alg.parS;parS.nErode=4;parS.nDilate=8;parS.conComp=2;        
        rec.M=refineMask(mean(abs(rec.xS),4:ND),parS,rec.Geom.S.MS);
        rec.xS=gather(rec.xS);
    end
       
    if isfield(rec,'PS');MTM=rec.Geom.PS.MT;else MTM=rec.Geom.S.MT;end
    if ~rec.Alg.useMasking;rec.Mx=[];    
    else rec.Mx=mapVolume(rec.M,rec.x,MTM,rec.Geom.x.MT,[],[],0,'nearest');
    end;rec.M=gather(rec.M);

    NF=rec.Enc.x.AcqN(2);
    F=build1DFTM(NF);
    rec.F=cell(1,NP);rec.G=cell(1,NP);
    for n=1:NP
        rec.F{n}=F(:,rec.Ay==n);
        if ~isempty(rec.Az);rec.G{n}=dynInd(rec.Az,find(rec.Ay==n),2);end
    end
    m=mod(NF,2);
    rec.y=encodeMSAlignedSense(((-1)^m)*rec.x,[],rec.F);
    if ~isempty(rec.Az)
        for n=1:NP;rec.y{n}=dynInd(rec.y{n},1:size(rec.y{n},3)/rec.Enc.(t).MultiBandFactor,3);end
    end
    rec.x=gather(rec.x);
    F=build1DFTM(NF);
    F=fftGPU(F,1);
    F=shifting(F,{m});
    F=ifftGPU(F,1);
    for n=1:NP;rec.F{n}=F(:,rec.Ay==n);end
    rec.x=conjugateGradientMSAlignedSense(rec.y,rec.Sx,rec.F,rec.G,rec.Mx);
    if ~isempty(rec.Alg.gibbsRinging)
        if any(rec.Alg.gibbsRinging>0);rec.x=gibbsRingingFilter(abs(rec.x),1,rec.Alg.gibbsRinging).*signz(rec.x);end%Appears safer to operate in magnitude for gradient echo
        rec.x=gather(rec.x);
    else
        rec.x=gather(rec.x);
        rec.x=unring(abs(rec.x)).*signz(rec.x);
    end

    for n=1:NP;rec.y{n}=gather(rec.y{n});end
    if isfield(rec,'Sx');rec=rmfield(rec,'Sx');end
    if isfield(rec,'Mx');rec=rmfield(rec,'Mx');end
    if isfield(rec,'PSx');rec=rmfield(rec,'PSx');end

%%%HEREHERHERE---GATHERSTRUCT AND SAVE (AFTER DELETING X)

    %return 
%     %rec.x=sum(conj(rec.S).*rec.x,4)./(sum(abs(rec.S).^2,4)+1e-12);
%     %return
%     
%     %RECONSTRUCT DISORDER
%     
%     %PARSE TO PHILIPS
%     %NAMES
%     rec.Names.Name=rec.Nam.dataIn;
%     rec.Names.pathOu=rec.Nam.caseIn;
%     rec.Plan.Suff='';
%     rec.Plan.SuffOu='';
%     %DATA
%     rec.y=rec.x;
%     TW=rec.TW.x;
%     if useGPU
%         rec.y=gpuArray(rec.y);rec.S=gpuArray(rec.S);rec.M=gpuArray(rec.M);
%     end
%     rec.S=permute(rec.S,[1 3 2 4]);rec.y=permute(rec.y,[1 3 2 4]);rec.M=permute(rec.M,[1 3 2 4]);    
%     %GEOMETRY (MAY NEED REVISION)
%     rec.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta([1 3 2]);
%     rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);    
%     rec.Par.Mine.MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];
%     %rec.Par.Mine.APhiRec=rec.Geom.x.MT*inv(rec.Geom.x.Co);
%     rec.Par.Mine.APhiRec=inv(rec.Par.Mine.MTT)*rec.Geom.S.MT(:,[1 3 2 4]);    
%     [~,rec.Par.Mine.Arot]=factorizeHomogeneousMatrix(rec.Par.Mine.APhiRec);
%     %TRAJECTORIES
%     rec.Par.Labels.TFEfactor=length(TW.image.Par);%This should be refined for shot-based sequences   
%     if rec.TW.x.hdr.Config.TurboFactor>1%We assume it is not steady state but based on shots
%         rec.Assign.TimeStamp=rec.TW.x.image.timestamp-rec.TW.x.image.timestamp(1);       
%         meduntimestamp=median(unique(diff(rec.Assign.TimeStamp)));%We look for timestamps larger than the median to identify boundaries between shots
%         shotsizes=diff([0 find(diff(rec.Assign.TimeStamp)>meduntimestamp) length(rec.Assign.TimeStamp)]);       
%         rec.Assign.Shots=repelem(1:length(shotsizes),shotsizes);
%     end
%     rec.Par.Labels.ZReconLength=1;            
%     NY=size(rec.y);NY(end+1:4)=1;
%     rec.Assign.z{2}=TW.image.Par-NY(2)/2-1;
%     rec.Assign.z{3}=permute(TW.image.Lin-NY(3)/2-1,[1 3 2]);            
%     %SEQUENCE PARAMETERS
%     rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;
%     rec.Par.Labels.TE=TW.hdr.MeasYaps.alTE{1}/1000;
%     rec.Par.Labels.FlipAngle=TW.hdr.MeasYaps.adFlipAngleDegree{1};
%     rec.Par.Labels.ScanDuration=TW.hdr.MeasYaps.lScanTimeSec;
%  
%     %DISORDER PARAMETERS
%     rec.Alg.WriteSnapshots=1;%To write snapshots
%     rec.Alg.parXT.groupSweeps=1;%Factor to group sweeps
%     rec.Alg.parXT.perc=[0.9 0.9 0.9];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
%     %rec.Alg.parXT.tolerSolve=1e-4;%Energy decrease in last iteration
%     %rec.Alg.parXT.saveFinal=1;%Save final results for further inspection
%     %rec.Alg.parXT.exploreMemory=1;%To explore convergence without running the main methods
%     
%     %RUN DISORDER
%     rec=solveXT(rec);
%     
%     %REVERT PERMUTATIONS
%     rec.x=permute(rec.x,[1 3 2 4]);
end

%PARTIAL FOURIER

%REMOVE OVERSAMPLING

%GEOMETRY CORRECTION
