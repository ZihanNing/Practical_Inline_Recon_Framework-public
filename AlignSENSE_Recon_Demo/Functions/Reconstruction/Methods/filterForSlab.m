

function [filterSize, isSlab, idx] = filterForSlab(FOVmm)


factor = .3;

isLowFOV = FOVmm<factor*multDimMax(FOVmm);
isSlab = any(isLowFOV);

if isSlab 
    if multDimSum(isLowFOV)==1
        warning('Slab in 2 direction is not expected and hence filter is disabled.');
        filterSize = [20 20 20];
    end
    
    idx = find(isLowFOV==1);
    filterSize = [20 20 20];
    filterSize(idx)=1;
else
    filterSize = [20 20 20];
end
