function grEn=gradientEntropy(x)

%GRADIENTENTROPY  Computes the gradient entropy of a multidimensional
%array
%   GREN=GRADIENTENTROPY(X)
%   * X is the multdimensional array
%   ** GREN is the gradient entropy
%

ND=numDims(x);
grEn=0;
for n=1:ND
    xd{n}=abs(x-dynInd(x,[2:size(x,n) 1],n));
    xd{n}=xd{n}+1e-12;
end
xd=cat(ND+1,xd{:});
xd=xd/(sum(xd(:)));    
grEn=grEn-gather(sum(xd(:).*log2(xd(:))))/log2(numel(xd));%Maximum possible entropy, all outcomes equally probable
