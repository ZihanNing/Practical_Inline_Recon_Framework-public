
function pGrouped = groupPT (p, stateSample, NStates)

% PT data as  NCha x NSamples
% stateSample as 1 x NSamples
% pGrouped as NCha x NStates

if nargin<3 || isempty(NStates);NStates = max(stateSample(:));end
if isa(stateSample, 'gpuArray'); stateSample=gather(stateSample);end

N = size(p);
NCha = N(1);
NSamples = N(2);

pGrouped = NaN([NCha, max(stateSample(:))],'like', p);%Fill with max(stateSample) since this will be the dimension of what comes out of regionprops
isComplex =  ~isreal(p);

for i = 1:NCha%each coil individually
    
    ttReal = regionprops(stateSample,dynInd(real(p),i,1),'MeanIntensity');ttReal={ttReal.MeanIntensity};ttReal=cat(2,ttReal{:});

    if isComplex
        ttImag = regionprops(stateSample,dynInd(imag(p),i,1),'MeanIntensity');ttImag={ttImag.MeanIntensity};ttImag=cat(2,ttImag{:});
        pGrouped = dynInd(pGrouped,i,1,ttReal + 1i* ttImag); 
    else
        pGrouped = dynInd(pGrouped,i,1,ttReal); 
    end
    
end

%Set empty shots to zero PT signal
if isComplex; valToFill = NaN + 1i*NaN; else; valToFill=NaN;end
pGrouped = dynInd(pGrouped,size(pGrouped,2)+1:NStates,2,valToFill);%Extend with NaNs at the end where you might have zero-samples shots
idxEmpty = isnan(pGrouped(1,:));

fillType=0;
if fillType==1%Fill with zero
    pGrouped = dynInd(pGrouped,idxEmpty,2,0);
elseif fillType==2%Remove NaNs
    pGrouped = dynInd(pGrouped,~idxEmpty,2);
end


