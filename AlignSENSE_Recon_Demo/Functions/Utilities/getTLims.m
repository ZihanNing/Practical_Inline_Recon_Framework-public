

function [limPlot] = getTLims(T, MS, addError)

if nargin<2 || isempty(MS);MS=1;end
if nargin<3 || isempty(addError);addError=1;end

limPlot = 1.1*[max(MS)*repmat(multDimMax(abs(dynInd(T,1:3,6))),[1 2]) ;
               repmat(convertRotation(multDimMax(abs(dynInd(T,4:6,6))),'rad','deg'),[1 2]) ] .* repmat([-1 1],[2 1]);

if addError;limPlot=limPlot+repmat([-1e-3 1e-3],[2,1]);end

limPlot = gather(limPlot);