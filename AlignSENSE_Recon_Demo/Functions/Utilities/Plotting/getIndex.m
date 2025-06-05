

function [idxPlotND] = getIndex(N,idx,supFOV,MT,isRASId)

%%% Extract FOV and do NIFTI index handling
for i=1:3; if isempty(supFOV{i}); supFOV{i} = 1:N(i);end;end
if isempty(idx); idx = ceil((N(1:3)+1)/2);end
idxPlotND = idx -  ([supFOV{1}(1) supFOV{2}(1) supFOV{3}(1)]-1);

if ~isempty(MT)
    [perm, fl] = T2perm(MT(1:3,1:3));%Flips and permutes applied to the array to make logical axes approximate the RAS axes.

    if isRASId==0
        idPlot(fl==1) = (N(fl==1)+1) - idPlot(fl==1);    
        idPlot = idPlot(perm(1:3));
    end
end