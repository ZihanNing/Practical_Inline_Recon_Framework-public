function visResidualsOld(ry,pau)

%VISRESIDUALS   Visualizes reconstruction residuals
%   VISRESIDUALS(RY)
%   * RY are the residuals to be visualized
%

if nargin<2 || isempty(pau);pau=1;end

FontSizeA=30;
figure
subtightplot(2,2,1,[0 0.075])
imshow(log(fftshift(multDimSum(ry,3:5))),[])
title('Log-error (spectrum)','FontSize',FontSizeA)
subtightplot(2,2,2,[0 0.075])
A=ry~=0;
plot(1:size(ry,5),squeeze(multDimSum(ry,1:3)./multDimSum(A,1:3))','*')
ylabel('Mean residual','FontSize',FontSizeA);
xlabel('Motion state','FontSize',FontSizeA);

[~,dV]=estimatePowerLaw(sqrt(ry),1,1.9);
subtightplot(2,2,4,[0 0.075])
plot(1:size(ry,5),-dV(:)','*')
ylabel('Power law','FontSize',FontSizeA);
xlabel('Motion state','FontSize',FontSizeA);

[~,dV,xn,Wn,p]=estimatePowerLaw(sqrt(sum(ry,5)),1,1.9);
subtightplot(2,2,3,[0 0.075])
plot(Wn(:),xn(:),'*')
hold on
plot(Wn(:),p(1)+Wn(:)*p(2),'*');
ylabel('Log-energy','FontSize',FontSizeA);
xlabel('Log-discrete-grid','FontSize',FontSizeA);
fprintf('Estimated law: %.2f\n',dV);

set(gcf, 'Position', get(0,'Screensize'))  

if pau;pause;end