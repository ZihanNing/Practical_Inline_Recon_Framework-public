function useGPUFlag=useGPU

%USEGPU   Allows to block GPU computations by setting the corresponding 
%flag to 0
%   USEGPUFLAG=USEGPU
%   ** USEGPUFLAG is a flag to perform GPU computations
%

useGPUFlag=1;
useGPUFlag=single(useGPUFlag && gpuDeviceCount);
