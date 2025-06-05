function visGhost(P)

%VISGHOST   Visualizes Nyquist ghosting correction information
%   VISGHOST(P)
%   * P is the Nyquist ghost correction data
%

N=size(P);N(end+1:3)=1;
P=reshape(P,[N(1:2) prod(N(3:end))]);

FontSizeA=30;
figure
for n=1:2
    Ps=dynInd(P,n,2);
    Ps=permute(Ps,[1 3 2]);
    subtightplot(1,4,2*n-1,[0 0])
    imshow(angle(Ps),[-pi pi])
    ylabel('Readout','FontSize',FontSizeA);
    xlabel('Phase','FontSize',FontSizeA)
    subtightplot(1,4,2*n,[0 0])
    imshow(abs(Ps),[])
    ylabel('Readout','FontSize',FontSizeA);
    xlabel('Magnitude','FontSize',FontSizeA)    
end
set(gcf, 'Position', get(0,'Screensize'))  
pause