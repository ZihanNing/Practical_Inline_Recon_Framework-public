

function E = convertChiToTaylorMaps (E)

if isfield(E,'Ds')
    
    if ~isfield(E.Ds,'chi'); warning('convertChiToTaylorMaps:: aborted since no susceptibility distribution chi provided.');return;end
    
    %%% FOURIER DOMAIN
    CHI = E.Ds.chi;%Already padded
    %CHI = padArrayND(CHI, E.Ds.padDim,1,0,'both');
    for i=1:3; CHI=fftGPU(CHI, i);end

    %%% CREATE DIFFERENTIAL DIPOLE KERNELS
    H0 = E.Ds.H0;%In T
    gamma = E.Ds.gamma;%MHz/T

    kernels = [];%Already padded
    if ~isempty(E.Ds.kernelStruct.Kequilibrium); kernels = cat(6, kernels, E.Ds.kernelStruct.Kequilibrium);end
    if ~isempty(E.Ds.kernelStruct.Kroll); kernels = cat(6, kernels, E.Ds.kernelStruct.Kroll);end
    if ~isempty(E.Ds.kernelStruct.Kpitch); kernels = cat(6, kernels, E.Ds.kernelStruct.Kpitch);end

    kernels = H0 * gamma * kernels;%Hz/rad for LCMs and in Hz for B0
    
    %%% FILTER KERNEL
    kernels = bsxfun(@times, fftshift(E.Ds.kernelStruct.H), kernels);%Filter kernel to avoid gibs artefacts
    for i=1:3; kernels = ifftshiftGPU(kernels,i);end %CHI not shifted so shift back

    %%% APPLY KERNEL
    kernels = bsxfun(@times, kernels, CHI);

    %%% MOVE KERNELS TO IMAGE DOMAIN
    for i=1:3; kernels=ifftGPU(kernels, i); end
    kernels = real(kernels);
    
    %INVERSE PADDING
    kernels = padArrayND(kernels,E.Ds.padDim,0,[],'both');

    %%%EXTRACT LCMs AND B0 FIELD
    if ismember(0,E.Ds.orderSusc)
        E.Ds.B0 = dynInd(kernels,1,6);%B0
        if ismember(1,E.Ds.orderSusc)
            E.Ds.D = dynInd(kernels,2:3,6);%LCMs
        end
    elseif ismember(1,E.Ds.orderSusc)
        E.Ds.D = dynInd(kernels,1:2,6);%LCMs
    end

else
   warning('convertChiToTaylorMaps:: No E.Ds structure!') 
end
