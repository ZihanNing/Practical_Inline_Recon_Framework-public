function x=LADODF(A,y,p,demed)

%LADODF  Estimates the ODFs up to a given order by nnIRWLS
%   X=LADODF(A,Y,{P},{DEMED})
%   * A is the DWI data
%   * Y are the corresponding diffusion measures
%   * {P} is the fitting norm
%   * {DEMED} balances for possible skewness in the data
%

if ~exist('demed','var') || isempty(demed);demed=0;end


if ~exist('p','var') || isempty(p);p=2;end
if ~exist('nIt','var') || isempty(nIt);nIt=300;end
if ~exist('tol','var') || isempty(tol);tol=1e-3;end
if ~exist('de','var') || isempty(de);de=1e-3;end%1e-6 is unstable

if p==2
    x=lsqnonneg(A,y);
else
    Eprev=inf;    
    res=y;res(:)=1;
    pon=2;%Initial homotopy value
    redN=0.8;%Reduction factor to slowly approach the target norm
    for n=1:nIt
        x=lsqnonneg(bsxfun(@times,res,A),res.*y);
        res=A*x-y;     
        E=sum(abs(res).^p);
        %if demed;res=res+median(res);end
        res=abs(res);
        if (Eprev-E)/Eprev<tol && n>1           
            break;
        else
            pon=max(p,redN*pon);
            res=(res+de).^(pon-2)/2;
            Eprev=E;                        
        end          
    end
    if n==nIt;fprintf('LAD reached the maximum number of iterations\n');end
end

%function histRes(r)
%[Fi,xi]=ksdensity(r);
%figure            
%plot(xi,Fi)
%hold on
%plot(median(r),0.01,'r*')
%plot(mean(r),0.01,'g*')
%grid on
%end
%