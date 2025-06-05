function [c,a,W]=robustLossEstimation(x,w)

%ROBUSTNESSESTIMATION   Estimates the robust function to be applied for
%data X. It uses the formulation in [1] JT Barron, "A General and Adaptive 
%Robust Loss Function," ArXiv:1701.03077v6:1-19.
%   [C,A,W]=ROBUSTNESSESTIMATION(X,{W})
%   * X are the input samples
%   * {W} are a set of reliability weights for the input samples. It
%   defaults to 1
%   * C is the scale parameter for robustness
%   * A is the shape parameter for robustness
%   * W are the weights to be applied in iteratively reweighted least
%   squares
%

if nargin<2 || isempty(w);w=ones(size(x),'like',x);end

c=sqrt(sum((x(:).^2).*w(:))./sum(w(:)));
a=2;

e=1e-5;%Tolerance factor

for n=1:20    
    a=fminbnd(@(af)negativeLogLikelihoodRobustLoss(x,af,c,w,e),0,2);
    c=fzero(@(cf)negativeLogLikelihoodCDerivativeRobustLoss(cf,x,a,w,e),c);
end

if nargout>=3;W=single(weightsRobustLoss(x,a,c,e));end


function [y,b]=negativeLogLikelihoodCDerivativeRobustLoss(c,x,a,w,e)
    if nargin<4 || isempty(w);w=ones(size(x),'like',x);end
    if nargin<5 || isempty(e);e=1e-5;end

    b=abs(2-a)+e;
    d=a+e-2*e*(a<0);
    xc2=(x/c).^2;
    xc2b=xc2/b;
    xc2b1=xc2b+1;
    y=(1/c)*(1-(xc2b1.^(d/2-1)).*xc2);
    y=sum(w(:).*y(:))./sum(w(:));
end


function y=negativeLogLikelihoodRobustLoss(x,a,c,w,e)
    if nargin<4 || isempty(w);w=ones(size(x),'like',x);end
    if nargin<5 || isempty(e);e=1e-5;end

    y=robustLoss(x,a,c,e)+log(c)+logZ1(a);
    y=sum(w(:).*y(:))./sum(w(:));

    function b=logZ1(a)
        xa=4*a;
        ps=[1.4913 1.36416 1.29171 1.23529 1.1855 1.13753 1.08719 1.02749 0.91894];
        ms=[-0.19877 -0.08884 -0.0633 -0.05293 -0.04901 -0.04958 -0.0556 -0.07346 -0.20015];
        i0=floor(min(max(xa,0),length(ps)-2))+1;
        p=[ps(i0:i0+1) ms(i0:i0+1)];
        t=xa-(i0-1);
        h(2)=t.*t.*(-2*t+3);
        h(1)=1-h(2);
        h(4)=t.*t.*(t-1);
        h(3)=h(4)+t.*(1-t);

        if t<0
            b=ms(1)*t+ps(1);
        elseif t>1
            b=ms(end)*(t-1)+ps(end);
        else
            b=sum(p.*h);
        end
    end
    
    function y=robustLoss(x,a,c,e)
        if nargin<4 || isempty(e);e=1e-5;end%Tolerance I guess
        b=abs(2-a)+e;
        d=a+e-2*e*(a<0);
        y=(b/d)*((((x/c).^2)/b+1).^(d/2)-1);    
    end
end

function y=weightsRobustLoss(x,a,c,e)
    if nargin<4 || isempty(e);e=1e-5;end    
    
    b=abs(2-a)+e;
    d=a+e-2*e*(a<0);
    y=(c^(-2))*(((x/c).^2/b+1)).^(d/2-1);
end

end

% %DERIVATIVE OF NEGATIVE LOG LIKELIHOOD WITH RESPECT TO A, NOT TESTED FOR COMPLIANCE!
% function y=daGeneralRobustLoss(x,a,c,e)
% 
% if nargin<4 || isempty(e);e=1e-5;end%Tolerance I guess
% 
% %Function to compute the general robust loss based on A General and 
% %Adaptive Robust Loss Function, Jonathan T. Barron, Google
% 
% b=abs(2-a)+e;
% d=a+e-2*e*(a<0);
% 
% xc2=(x/c).^2;
% xc2b=xc2/b;
% xc2b1=xc2b+1;
% grl=robustLoss(x,a,c,e);
% 
% y=(1/b)*(grl-0.5*xc2.*(xc2b1.^(d/2-1)))*sign(a-2)+(1/d)*(0.5*b*log(xc2b1).*(xc2b1.^(d/2))-grl);
% if a>2-e;d1=2-e;d2=2;
% elseif a<e;d1=0;d2=e;
% else d1=a-e/2;d2=a+e/2;
% end
% [~,b2]=negativeLogLikelihoodRobustLoss(x,d2,c,e);
% [~,b1]=negativeLogLikelihoodRobustLoss(x,d1,c,e);
% y=y+(b2-b1)/e;

% %DERIVATIVE OF NEGATIVE LOG LIKELIHOOD WITH RESPECT TO C, ALREADY INTRODUCED!
% function y=dcGeneralRobustLoss(x,a,c,e)
% 
% if nargin<4 || isempty(e);e=1e-5;end%Tolerance I guess
% 
% %Function to compute the general robust loss based on A General and 
% %Adaptive Robust Loss Function, Jonathan T. Barron, Google
% %This is the derivative of the negative log likelihood
% 
% b=abs(2-a)+e;
% d=a+e-2*e*(a<0);
% 
% xc2=(x/c).^2;
% xc2b=xc2/b;
% xc2b1=xc2b+1;
% 
% y=(1/c)*(1-(xc2b1.^(d/2-1)).*xc2);

