function [y,z,hd]=phaseDistribution(x,K,typ)

%PHASEDISTRIBUTION returns the phase distribution of complex additive 
%Gaussian noise
%   [Y,Z]=PHASEDISTRIBUTION({X},{K},{TYP}) returns the phase distribution and moments of
%   complex additive Gaussian noise
%   X is the array where to compute the distribution (rows)
%   K is the Rice factor (s^2/(2*sigma^2)) (column) 
%   Defaults to 1.
%   TYP are the types of moments to compute, 0 for trigonometric moments
%   (cosine mean / (sine variance+cosine variance)/2) (default), 1 for
%   angular moments (variance)
%   Y is the corresponding distribution
%   Z are the default moments   
%   HD is the Hellinger distance with regard to the uniform phase 
%   distribution
%

if nargin<1 || isempty(x);x=-pi:0.01:pi;end
if nargin<2 || isempty(K);K=1;end
if nargin<3 || isempty(typ);typ=0;end

x=x(:);
K=K(:)';
K(K>1e8)=1e8;
Ksqco=bsxfun(@times,sqrt(K),cos(x));

%y=bsxfun(@times,(exp(-K)/(2*pi)),1+sqrt(pi)*Ksqco.*exp(Ksqco.^2).*erfc(-Ksqco));
%We need to insert the exponential inside for numerical stability:
y=bsxfun(@plus,exp(-K),sqrt(pi)*Ksqco.*exp(bsxfun(@minus,Ksqco.^2,K)).*erfc(-Ksqco))/(2*pi);

if nargout>1 
    x=-pi:0.01:pi;x=x(:);
    Ksqco=bsxfun(@times,sqrt(K),cos(x));
    if typ==0
        z{1}=cos(x);%Mean of real noise component
        z{2}=cos(x).^2;%Second moment of real noise component
        z{3}=sin(x).^2;%Second moment of imaginary noise component
        for n=1:length(z)
            z{n}=bsxfun(@times,z{n},bsxfun(@plus,exp(-K),sqrt(pi)*Ksqco.*exp(bsxfun(@minus,Ksqco.^2,K)).*erfc(-Ksqco))/(2*pi));
            z{n}=sum(bsxfun(@times,z{n},gradient(x)),1);
        end
        z{2}=z{2}-z{1}.^2;%Variance
        z{2}=(z{2}+z{3})/2;%Averaged variance
        z(3)=[];
    else
        z{1}=x.^2;
        for n=1:length(z)
            z{n}=bsxfun(@times,z{n},bsxfun(@plus,exp(-K),sqrt(pi)*Ksqco.*exp(bsxfun(@minus,Ksqco.^2,K)).*erfc(-Ksqco))/(2*pi));
            z{n}=sum(bsxfun(@times,z{n},gradient(x)),1);
        end
    end
end
if nargout>2
    x=-pi:0.01:pi;x=x(:);
    Ksqco=bsxfun(@times,sqrt(K),cos(x));
    hd=sqrt(bsxfun(@plus,exp(-K),sqrt(pi)*Ksqco.*exp(bsxfun(@minus,Ksqco.^2,K)).*erfc(-Ksqco)))/(2*pi);
    %hd=sqrt(max(1-sum(bsxfun(@times,hd,gradient(x)),1),0));
    hd=1-sum(bsxfun(@times,hd,gradient(x)),1);
    if K(1)==0;hd=sqrt(hd-hd(1));else hd=sqrt(max(1-sum(bsxfun(@times,hd,gradient(x)),1),0));end
end
    