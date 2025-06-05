

function TRes = resampleT(T, NStatesNew, interpolationMethod)
%SEE UPSAMPLET.m from Luclio's repo
N = size(T);
nD = ndims(T);

if nargin < 2 || isempty(NStatesNew); NStatesNew=size(T,nD-1);end
if nargin < 3 || isempty(interpolationMethod); interpolationMethod='cubic';end

NStates = size(T,nD-1);
NMotParam = size(T,nD);
NNew = N; NNew(nD-1)=NStatesNew;
TRes = zeros(NNew,'like',T);

x = linspace(0,1,NStates);%0 and 1 don't really matter - to be more accurate, use the actual time stamps in x
xQuery = linspace(0,1,NStatesNew);

for i=1:NMotParam
    tt = interp1(x,squeeze(dynInd(T,i,6)),xQuery,interpolationMethod);
    TRes = dynInd( TRes, i,6, tt);
end