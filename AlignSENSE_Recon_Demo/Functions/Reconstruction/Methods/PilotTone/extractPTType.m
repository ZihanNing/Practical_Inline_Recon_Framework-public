

function p = extractPTType(p, signalUsage )

isSlice = size(p,4)>1;
if isSlice; dimCha=4;dimSamples=1:2; else; dimCha=1;dimSamples=2;end

%%% Class of signal
if signalUsage.useMagn==1
    p = abs(p); 
elseif signalUsage.usePhase==1
    p= angle(p); 
elseif signalUsage.useRealImag == 1
    p = cat(dimCha, real(p), imag(p));
end

