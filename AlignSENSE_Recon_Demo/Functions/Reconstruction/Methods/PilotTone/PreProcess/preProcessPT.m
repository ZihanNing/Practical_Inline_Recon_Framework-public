
function [PT, parPT] = preProcessPT(PT, parPT, NY, kIndex, deb)

if nargin<3 || isempty(NY);NY=[];end
if nargin<4 || isempty(kIndex);kIndex=[];end
if nargin<5 || isempty(deb);deb=1;end

%%% SET ASIDE VARIABLES
pSlice = PT.pSlice;
PT.NChaPTAcq= size(pSlice,4);% # Non-compressed channels
if parPT.signalUsage.useRealImag; PT.NChaPTRec = 2*size(pSlice,4);else; PT.NChaPTRec = size(pSlice,4) ;end%Will be over-written when you do svd analysis
reOrderTemporally = parPT.signalUsage.eigTh<1 || parPT.PTFilter.medFiltKernelWidthidx>1; 
reOrderTemporally = reOrderTemporally || ~isempty(parPT.signalUsage.orderPreProcessing) ;%only for debugging during development
fprintf('<strong>Pilot Tone</strong>: Pre-processing\n')

%%% MULTI-BAND HANDLING
if ~isempty(PT.idxMB)
    if parPT.signalUsage.combineMB==1%Extract peak in RO direction
        pSlice = dynInd(pSlice, PT.idxMB,3);
    elseif parPT.signalUsage.combineMB==2%Average frequencies
        pSlice = multDimMea(pSlice, 3);
    else
       assert(size(pSlice,3)==1, 'Reconstruction not able to work with multi-band Pilot Tone signal.');
    end
    fprintf('	Removing multi-band PT signal.\n')
end

%%% RE-ORDER IN TEMPORAL DOMAIN IF NEEDED
if reOrderTemporally
    %%% Move the PT signal to the hybrid space where samples have meaning
    for n=1:2
        %pSlice=fftshiftGPU(pSlice,n);%Already done in solveXTB_PT script.
        pSlice=fftGPU(pSlice,n)/(NY(n));
        pSlice=fftshiftGPU(pSlice,n);
    end
    %%% Order in time
    idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling
    pTime=[];
    for i=1:NY(5)%Multiple repeats
        idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
        pTime=cat(1,  pTime, ...
                      dynInd( resPop(dynInd(pSlice,i,5) ,1:2,[],1),...
                              dynInd(idx,idxToExtract,2) ,1));
    end
    %pTime = pSlice %This was old one when not accounting for multiple averages
    %pTime = resPop(pTime,1:2,[],1) %This was old one when not accounting for multiple averages
    %pTime = dynInd(pTime,idx,1) %This was old one when not accounting for multiple averages
else
    pTime = pSlice;
end
%if reOrderTemporally;fprintf('		Conditioning number = %d\n', cond(permute(pTime,[1 4 2:3]).'));end

if deb && reOrderTemporally
    tt=permute(pTime,[1 4 2:3]).';
    visPTSignal(tt,'pTime',[],[],2,[],[],89,'Original Pilot Tone signal');
end

%%% TIME-INDEPENDENT PRE-PROCESSING
for step = parPT.signalUsage.orderPreProcessing
    if step == 1 %Removing complex mean values along PE directions
        pTime = bsxfun(@minus, pTime, multDimMea((pTime),1:2) .* exp(0*1i*multDimMea(angle(pTime))) );
        fprintf('	Removing complex mean values along PE directions.\n')
    end
    if step == 2 %Removing mean phase along PE directions
        pTime = bsxfun(@rdivide, pTime, sign(multDimMea(pTime,1:2)));
        fprintf('	Removing mean phase along PE directions.\n')
    end
    if step == 3 %Normalising along coil dimension
        pTime = bsxfun(@rdivide, pTime , sqrt(normm(pTime,[],4)) );
        fprintf('	Normalising along coil dimension.\n')
    end
    if step == 4 %Whitening complex signal
        permWhit = [1 4 2:3];
        percSampes=[];
        pTime = ipermute(whiten(permute(pTime,permWhit),percSampes),permWhit);
        fprintf('	Whitening complex signal.\n')
    end
    if step == 5 %Making unit-variance signal
        pTime = bsxfun(@rdivide, pTime, std((pTime),[],4)+0*std(imag(pTime),[],4));
        fprintf('	Making unit-variance signal.\n')
    end
    %if reOrderTemporally;fprintf('		Conditioning number = %d\n', cond(permute(pTime,[1 4 2:3]).'));end
end

if deb && reOrderTemporally
    tt=permute(pTime,[1 4 2:3]).';
    visPTSignal(tt,'pTime',[],[],2,[],[],90,'After time-independent pre-processing');
    %[~,pReshaped] = reshapeP(rand([6 32]),permute(pTime,[1 4 2:3]).','mat2vec');
end
                
%%% DIMENSIONALITY REDUCTION USING SINGULAR VALUE DECOMPOSITION
if parPT.signalUsage.eigTh~=1
    if parPT.signalUsage.useRealImagSVD && ~parPT.signalUsage.useRealImag; warning('When using real/imag for SVD but not during calibration, there might be a channel filled with an unwanted component.');end
    if ~parPT.signalUsage.useRealImagSVD && ~parPT.signalUsage.useRealImag; warning('When using complex SVD but not during calibration, the channels might be mixed w.r.t to the singular values.');end
    
    %%% Perform SVD
    permSVD = [1 4 2:3];
    pTime = permute(pTime, permSVD);
    if parPT.signalUsage.useRealImagSVD; pTime = cat(2,real(pTime),imag(pTime)); end
    pTimeOrig = pTime;
    [U,S,V] = svd(pTime,'econ');
    
    %%% Select #singular vectors to keep
    sv = diag(S);
    if parPT.signalUsage.eigTh==0 %Automatic/tuned detection
        [~,nv]=screePoint(sv,.7);
    elseif parPT.signalUsage.eigTh < 0 && mod(parPT.signalUsage.eigTh,1)==0 %-eigTh is the pre-defined number of components
        nv = -parPT.signalUsage.eigTh;
    else %eigTh is fraction of sv to discard
        nv=find(sv>=sv(1)*parPT.signalUsage.eigTh,1,'last');
    end
        
    %%% Re-generate signal using set of singular values
    if parPT.signalUsage.useRealImagSVD
        nvSave = ceil(nv/2); %Needs to be even when we store it in complex number again
        pTime = U(:,1:2*nvSave);  
        pTime = pTime(:,1:nvSave) + 1i*pTime(:,nvSave+1:end);%Store as complex array for fft resampling functionality
        if parPT.signalUsage.useRealImag; PT.NChaPTRec = nv;else;PT.NChaPTRec = nvSave ;end
    else
        nvSave = nv;
        pTime = U(:,1:nvSave);
        if parPT.signalUsage.useRealImag; PT.NChaPTRec = 2*nv;else;PT.NChaPTRec = nvSave ;end
    end
        
    %%% Add a reporting stage with the number of componenents in each level

    pTime = ipermute(pTime, permSVD);
    
    if deb
        pTimeProj = U(:,1:nv)*S(1:nv,1:nv)*V(:,1:nv)';
        pTimeProjError = pTimeOrig-pTimeProj;
        visSVDThresholding(U,S,V,nv,91,[],0);
        visPTSignal(pTimeOrig.','pTime',[],[],[],[],[],93,'Original signal');
        visPTSignal(pTimeProj.','pTime',[],[],[],[],[],94,'Projected signal')        
        visPTSignal(pTimeProjError.','pTime',[],[],[],[],[],95,'Error')        
        %visPTSignal(pTimeProjError.'./(pTimeOrig.'+eps('single')),'pTime',[],[],[],[],[],96,'Relative Error')        
    end
    fprintf('	Low rank respresentation: using %d components.\n',nv);
    %fprintf('		Conditioning number = %d\n', cond(permute(pTime,[1 4 2:3]).'))
else
   if parPT.signalUsage.useRealImag; PT.NChaPTRec = 2*size(pTime,4);else;PT.NChaPTRec = size(pTime,4) ;end 
end

%%% FILTER IN TIME DOMAIN
medFiltKernelWidth = round(parPT.PTFilter.medFiltKernelWidthidx);
if medFiltKernelWidth>1
    gpu = isa(pTime,'gpuArray');
    if gpu; pTime = gather(pTime);end
    pTime = medfilt1( real(pTime),medFiltKernelWidth,[],1) + ...%real channel
            1i* medfilt1(imag(pTime),medFiltKernelWidth,[],1) ; %imaginary channel
    if gpu; pTime = gpuArray(pTime);end
    fprintf('	Median filtering with %-element kernel.\n',medFiltKernelWidth);
    fprintf('		Conditioning number = %d\n', cond(permute(pTime,[1 4 2:3]).'))
end    


%%% CONVERT TO SLICE DOMAIN
if reOrderTemporally
    %%% Create k-space array
    pSlice = zerosL(pSlice);
    pSlice = resPop(pSlice,1:2,[],1);%Flatten PE dimensions
    pSlice = dynInd(pSlice, 1:size(pTime,4),4);%Make coil dimensions compatible with pre-processing
    
    %%% Deal with multiple averages
    for i=1:NY(5)%Multiple repeats
        idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
        pSlice = dynInd(pSlice, {idx(idxToExtract), i} , [1 5] , dynInd(pTime, idxToExtract,1)  );
    end

    pSlice = resPop(pSlice,1,NY(1:2),1:2);%Reshape to PE dimensions
        
    %%% Move the PT signal back to the image domain
    for n=1:2
        pSlice=ifftshiftGPU(pSlice,n);
        pSlice=ifftGPU(pSlice,n)*(NY(n));
        %pSlice=fftshiftGPU(pSlice,n);
    end
else
    pSlice = pTime;
end

%%% ASSIGN NEW VARIABLES
PT.pSlice = pSlice;

