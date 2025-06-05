

function [ WePT, outlWePT, stateSample ] = getPTWeightingOld(pTime , mSt , NStates, th)

N = size(pTime);
NCha = N(1);
NSamples = N(2);

stateSample= permute( mSt(mSt<=NStates),[2 1]);

%take coil with highest signal
tt= multDimMea( abs(pTime),2);
[ ~, idxMax ] = max (tt,[],1);

p = dynInd(pTime,idxMax, 1);
p=abs(p);
isComplex =  ~isreal(p);

i=1; %Only one channel        
        
varReal = regionprops(stateSample,dynInd(real(p),i,1),'PixelValues');
varReal={varReal.PixelValues};varReal = cellfun(@(x) var(x(:)), varReal);

WePT = 1 - varReal./ max(varReal);
outlWePT = varReal>mean(varReal);

% if isComplex
%     ttImag = regionprops(stateSample,dynInd(imag(p),i,1),'PixelValues');
%     ttImag={ttImag.PixelValues};ttImag = cellfun(@(x) var(x(:)), ttImag);
%     pGrouped = dynInd(pGrouped,i,1,ttReal + 1i* ttImag); 
% else
%     pGrouped = dynInd(pGrouped,i,1,ttReal); 
% end
    