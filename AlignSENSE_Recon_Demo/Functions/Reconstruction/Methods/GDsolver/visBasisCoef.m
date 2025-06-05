
function [] = visBasisCoef(Db,n_plot,t, numFig)

if nargin<2 || isempty(n_plot); n_plot = 0:Db.SHorder;end
if nargin<3; t = [];end
if nargin<4; numFig = 777;end

assert(max(n_plot)<=Db.SHorder, 'Order of coefficients to plot should be smaller/equal to Db.SHorder');

NStates = size(Db.c,5);%E.NMs;
t = multDimMea(t,2:ndims(t));
if ~isempty(t);tplot =t; xlab='$t$ [s]';else tplot=1:NStates;xlab='$t$ [a.u]';end

FontSizeA=10;
ModLab=5;
LineWidth=3;
invCol=1;

h= figure(numFig);
set(h,'color','w')
for i=1:length(n_plot)
    n=n_plot(i);
    cidx = numCoefs(n-1)+(1:2*n+1); 

    c = dynInd(Db.c,cidx,6);% can also pick cr
    c = squeeze(c);

    subtightplot(1,length(n_plot),i,[0.02 0.02],[0.07 0.05],[0.03 0.05]) 
    plot(tplot, c,'LineWidth',LineWidth,'LineStyle','-')
    
    %set(gca,'FontSize',FontSizeA+ModLab)
    xlabel(xlab,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    ylabel('Amplitude [a.u.]','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSizeA+ModLab)
    title(sprintf('Solid Harmonics order %d',n),'Interpreter','latex','FontSize',FontSizeA+ModLab)
    
    
    labels = {};
    for m=-n:n; labels = cat(2,labels , sprintf('$SH_{%d,%d}$',n,m) );end
    leg = legend(labels)    ;
    set(leg,'Interpreter','latex','FontSize',FontSizeA+2);
    legend show
    axis tight

end


end


%%%FUNCTIONS
function [N]=  numCoefs(n)
N=0;
for i=0:n; N =N+2*i+1;end
end
