
function pSlice = preProcessPT(pSlice, parPT, idxMB, NY, kIndex)

%MB handling
if ~isempty(idxMB)
    if parPT.signalUsage.combineMB==1
        pSlice = dynInd(pSlice, idxMB,3);%Extract peak in RO direction
    elseif parPT.signalUsage.combineMB==2
        pSlice = multDimMea(pSlice, 3);%Average frequencies
    else
       assert( size(pSlice,3)==1, 'Reconstruction not able to work with multi-band Pilot Tone signal.');
    end
end

%Whitening phase  and Normalisation
if parPT.signalUsage.whiten || parPT.signalUsage.whitenPhase || parPT.signalUsage.normaliseCoils || parPT.signalUsage.eigTh<1
    
    %%% Move the PT signal to the hybrid space where samples have meaning
    for n=1:2
        %pSlice=fftshiftGPU(pSlice,n);%Already done in solveXTB_PT script.
        pSlice=fftGPU(pSlice,n)/(NY(n));
        pSlice=fftshiftGPU(pSlice,n);
    end
    
    %%% Order in time
    idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling

    pTime = pSlice;
    pTime = resPop(pTime ,1:2,[],1);
    pTime = dynInd(pTime, idx, 1);
        
    %%% Remove mean phase from the samples in the hybrid space    
    for step = parPT.signalUsage.orderPreProcessing
        if step ==1 
            if parPT.signalUsage.whiten; pTime = bsxfun(@minus, pTime, multDimMea((pTime),1:2) .* exp(0*1i*multDimMea(angle(pTime))) );end
        end
        if step ==2 
            if parPT.signalUsage.whitenPhase; pTime = bsxfun(@rdivide, pTime, sign(multDimMea(pTime,1:2)));end
        end
        if step ==3 
            if parPT.signalUsage.normaliseCoils; pTime = bsxfun(@rdivide, pTime , sqrt(normm(pTime,[],4)) );end
        end
    end
    
    %%% Use lower order singular vectors to reduce dimensionality
    if parPT.signalUsage.eigTh<1
        %Perform SVD
        permSVD = [1 4 2:3];
        pTime = permute(pTime, permSVD);
        if parPT.signalUsage.useRealImagSVD; pTime = cat(2,real(pTime),imag(pTime)); end
        pTimeOrig = pTime;
        [U,S,V] = svd(pTime,'econ');
        D = diag(S);
        %Select #singular vectors to keep
        if parPT.signalUsage.eigTh==0 %Automatic/tuned detection
            [~,nv]=screePoint(D,.7);
        else %parPT.signalUsage.eigTh is fraction of sv to keep
            nv=find(D>=D(1)*parPT.signalUsage.eigTh,1,'last');
        end
                %%% PLOTTING PARAMETERS
                LineWidth = 2;
                FontSize = 20;

                %%% PLOT
                h=figure(90);
                set(h,'color','w');
                plot([1:length(D)],D,'LineWidth',2); hold on 
                plot([nv,nv],[0,D(1)],'r--','LineWidth',LineWidth)
                hLeg = legend({'Signular values','Threshold'}); set(hLeg,'Interpreter','latex','FontSize',FontSize);
                xlabel('Singular value [$\#$]','interpreter', 'latex','FontSize',FontSize);
                ylabel('Value [a.u.]','interpreter', 'latex','FontSize',FontSize);
                title('Singular values', 'interpreter', 'latex', 'FontSize',FontSize)
                axis tight

                h=figure(91);set(h,'color','w');
                np = ceil(sqrt(size(U,2)));
                for i=1:size(U,2); subplot(np,np,i); plot(1:length(U(:,i)),real(U(:,i).'));end

        %Re-generate signal using set of singular values
        if parPT.signalUsage.useRealImagSVD; nv = nv+mod(nv,2); end%Needs to be even when we store it in complex number again
        pTime = U(:,1:nv);
        if parPT.signalUsage.useRealImagSVD; pTime = pTime(:,1:nv/2) + 1i*pTime(:,nv/2+1:end);end%Store as complex array for fft resampling functionality
        
        pTimeProj = U(:,1:nv)*S(1:nv,1:nv)*V(:,1:nv)';
        pTimeProjError = pTimeOrig-pTimeProj;
        
plotPTSignals(pTimeOrig.',1);
plotPTSignals(pTimeProj.',2);
plotPTSignals(pTimeProjError.',3);

        pTime = ipermute(pTime, permSVD);
    end

    %%% Filter
    medFiltKernelWidth = round(parPT.PTFilter.medFiltKernelWidthidx);
    if medFiltKernelWidth>1
        gpu = isa(pTime,'gpuArray');
        if gpu; pTime = gather(pTime);end
        pTime = medfilt1( real(pTime),medFiltKernelWidth,[],1) + ... %real channel
                1i* medfilt1(imag(pTime),medFiltKernelWidth,[],1) ; %imaginary channel
        if gpu; pTime = gpuArray(pTime);end
    end    
    
    %%% Create array
    pSlice = zerosL(pSlice);pSlice = resPop(pSlice,1:2,[],1);pSlice = dynInd(pSlice, 1:size(pTime,4),4);
    pSlice = dynInd(pSlice, idx, 1, dynInd(pTime, 1:size(pTime,1),1)  );
    pSlice = resPop(pSlice ,1,NY(1:2),1:2);
        
    %%% Move the PT signal back to the image domain
    for n=1:2
        pSlice=ifftshiftGPU(pSlice,n);
        pSlice=ifftGPU(pSlice,n)*(NY(n));
        %pSlice=fftshiftGPU(pSlice,n);
    end
end

