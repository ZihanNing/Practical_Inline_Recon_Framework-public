function rec=reconSENSE_gadg(rec)
%     if useGPU;rec.S=gpuArray(rec.S);rec.Ay=gpuArray(rec.Ay);rec.Az=gpuArray(rec.Az);rec.x=gpuArray(rec.x);rec.N.x=gpuArray(rec.N.x);end
    useGPU = 1;
    ND=16;
    if useGPU;rec.S=gpuArray(rec.S);rec.Ay=gpuArray(rec.Ay);rec.Az=gpuArray(rec.Az);rec.x=gpuArray(rec.x);end % ZN: ignored noise (temporarily for gadgetron version)
%     % ZN: for gadg verison, noise pre-whitening is skipped
%     if ~isempty(rec.N.x)
%         [rec.x,n]=standardizeCoils(rec.x,rec.N.x);
%         rec.S=standardizeCoils(rec.S,n);
%     end
    NP=gather(max(rec.Ay(:)));%Number of polarities
    
    % interpolate low-res refscan to high-res as the high-res dataset
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
end