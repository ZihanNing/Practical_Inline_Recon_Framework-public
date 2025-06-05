
function TGrouped = groupT (T, stateSample,NStates)

if nargin<3 || isempty(NStates);NStates = max(stateSample(:));end
if isa(stateSample, 'gpuArray'); stateSample=gather(stateSample);end

N = size(T);
TGrouped = NaN([N(1:4), max(stateSample(:)) , N(6) ]);%Fill with max(stateSample) since this will be the dimension of what comes out of regionprops

for i = 1:N(6)
    tt = regionprops(permute(stateSample,[1 3:5 2 6]),dynInd(T,i,6),'MeanIntensity');tt={tt.MeanIntensity};tt=cat(5,tt{:});
    TGrouped = dynInd(TGrouped,i,6,tt);
end
    
% %Set empty shots to zero motion parameters
TGrouped = dynInd(TGrouped,size(TGrouped,5)+1:NStates,5,NaN);%Extend with NaNs at the end where you might have zero-samples shots
idxEmpty = isnan( dynInd(TGrouped,1,6) );

fillType=0;
if fillType==1%Fill with zero
    TGrouped = dynInd(TGrouped,idxEmpty,5,0);
elseif fillType==2%Remove NaNs
    TGrouped = dynInd(TGrouped,~idxEmpty,5);
end


