

function p = filterPT(p, filterStruct, NY, kIndex, type)

if nargin < 2 || isempty(filterStruct); filterStruct=[];end
if nargin < 3 || isempty(NY); NY=[];end
if nargin < 4 || isempty(kIndex); kIndex=[];end
if nargin < 5 || isempty(type); type='pSlice';end

gpu = isa(p, 'gpuArray');

if strcmp(filterStruct.Type,'median') && round(filterStruct.medFiltKernelWidthidx)>1
    medFiltKernelWidth = round(filterStruct.medFiltKernelWidthidx);
    
    if strcmp(type,'pSlice')
        pSlice = p;

        %%% Move the PT signal to the fourier domain
        for n=1:2
            %pSlice=fftshiftGPU(pSlice,n);
            pSlice=fftGPU(pSlice,n)/sqrt(NY(n));
            pSlice=fftshiftGPU(pSlice,n);
        end

        %%% Extract temporal signal using sampling patters
        idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling
    %       
    %                 timeMat = zeros([prod(NY(1:2)) 1]);
    %                 timeMat(idx) = 1:length(idx);
    %                 timeMat = resSub(timeMat, 1, NY(1:2)); 
    %                 hitMat = single(timeMat>0);

        pTime = pSlice;
        pTime = resPop(pTime ,1:2,[],1);
        pTime = dynInd(pTime, idx, 1);

        %%% Filter using median filter
        if gpu; pTime = gather(pTime);end
        pTime = medfilt1( real(pTime),medFiltKernelWidth,[],1) + ... %real channel
                1i* medfilt1(imag(pTime),medFiltKernelWidth,[],1) ; %imaginary channel

        if gpu; pTime = gpuArray(pTime);end

        %%% Create array
        pTimeNew = zerosL(pSlice);pTimeNew = resPop(pTimeNew ,1:2,[],1);
        pTimeNew = dynInd(pTimeNew, idx, 1, dynInd(pTime, 1:size(pTime,1),1)  );
        pSlice = resPop(pTimeNew ,1,NY(1:2),1:2);

        %%% Move the PT signal back to the image domain
        for n=1:2
            pSlice=ifftshiftGPU(pSlice,n);
            pSlice=ifftGPU(pSlice,n)*sqrt(NY(n));
            %pSlice=fftshiftGPU(pSlice,n);
        end
        p = pSlice;
        
    elseif strcmp(type,'pTime')
        pTime = p;
        %%% Filter using median filter
        if gpu; pTime = gather(pTime);end
        pTime = medfilt1( real(pTime),medFiltKernelWidth,[],2) + ... %real channel
                1i* medfilt1(imag(pTime),medFiltKernelWidth,[],2) ; %imaginary channel

        if gpu; pTime = gpuArray(pTime);end
        p=pTime;
    else
        error('filterPTSlice:: Filtering for type %s not implemented.',type)
    end
end