
function Re_seg = res_segments(Re, NMs, nEc)

    %Re = computeEnergy(y,x,E,[],EH,[],[],[],[],1)/numel(y);
    Re_seg = [];
    for s = 1:NMs; Re_seg = cat( 1, Re_seg, multDimSum( dynInd(Re,nEc(s)+1:nEc(s+1),1)));end
    %Re_seg2 = accumarray(E.mSt, Re);%might be an issue with segments/states offset (shutter)
    
end

