for l=1:3
    data = fftshiftGPU(data,l); 
    data = fftGPU(data,l);%*NY(l);  %YB: was ifft originally
    data = ifftshiftGPU(data,l);
end