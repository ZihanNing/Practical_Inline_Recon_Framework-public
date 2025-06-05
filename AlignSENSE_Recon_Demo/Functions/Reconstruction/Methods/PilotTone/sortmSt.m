

function [mStSorted] = sortmSt(mSt,NStates)

uniqueSates = sort(unique(mSt));
newStates = 1:length(uniqueSates);
if nargin<2 || isempty(NStates);NStates=length(uniqueSates);end
assert(NStates==length(uniqueSates),'sortmSt:: NStates specified does not agree with the provided mSt.')

mStSorted=zerosL(mSt);
for i=1:NStates
    idx = mSt==uniqueSates(i);
    mStSorted(idx)=newStates(i);
end
assert(~any(mStSorted==0),'sortmSt:: State 0 not allowed.')

% mStSorted=mSt;
% idxRunning =1;
% 
% for stateSorted=1:NStates
%    foundNext=0;
%    while foundNext==0
%        idx = find(mStSorted == stateSorted);
%        if isempty(idx)
%            mStSorted(idxRunning:end) = mStSorted(idxRunning:end)-1;
%        else
%            foundNext=1;
%            idxRunning = idxRunning+length(idx);
%        end
%    end
% end
% 
