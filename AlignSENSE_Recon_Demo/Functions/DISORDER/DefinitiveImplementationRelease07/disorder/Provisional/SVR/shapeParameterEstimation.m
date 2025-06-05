function [b,c]=shapeParameterEstimation(x,w)

%SHAPEPARAMETERESTIMATION   Estimates the shape parameter of a generalized
%Gaussian distribution using [1] KS Song, "A Globally Convergent and 
%Consistent Method for Estimating the Shape Parameter of a Generalized 
%Gaussian Distribution", IEEE Trans Inf Th, 52(2):510-527, 2006.
%   B=SHAPEPARAMETERESTIMATION(X,{W})
%   * X are the input samples
%   * {W} are a set of reliability weights for the input samples. It
%   defaults to all ones
%   * B is the shape parameter for robustness
%

if nargin<2 || isempty(w);w=ones(size(x),'like',x);end

b=2;
nMax=20;
tol=1e-4;
x=abs(x)+1e-9;
for n=1:nMax
    bold=b;
    m1=sum((x(:).^(1*b)).*w(:))./sum(w(:));
    m2=sum((x(:).^(2*b)).*w(:))./sum(w(:));
    m1log=sum((x(:).^(1*b)).*log(x(:)).*w(:))./sum(w(:));
    m2log=sum((x(:).^(2*b)).*log(x(:)).*w(:))./sum(w(:));
    bd=(m2/m1^2-(1+b));
    bdd=(2*(m2log*m1^2-m1log*m2*m1)-m1^4)/m1^4;
    b=b-bd/bdd;
    if abs(b-bold)<tol;break;end
end
c=(b*sum((x(:).^(1*b)).*w(:))./sum(w(:)))^(1/b);