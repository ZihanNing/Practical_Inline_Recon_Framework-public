function gpuInfo = checkGPUAvailability()
    gpuInfo = gpuDeviceCount;
    if gpuInfo > 0
        fprintf('There are %d GPUs available for computation.\n', gpuInfo);
    else
        error('No GPU available for computation.');
    end
end