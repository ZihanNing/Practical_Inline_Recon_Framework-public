function visSensitivities(S,visSenPar,voxsiz)

%VISSENSITIVITIES   Visualizes the compressed sensitivities at half of the
%readout direction, the compressed sensitivities at that location, and the
%(normalized) l2 and linf norm of the compressed decomposition of the 
%3D sensitivity maps 
%   VISSENSITIVITIES(S)
%   * S is the sensitivity map
%   * VISSENPAR indicates whether to write the data (2) or not (1)
%   * VOXSIZ are the voxel sizes of the coils. They default to all 1s
%

if nargin<2 || isempty(visSenPar);visSenPar=1;end
if nargin<3 || isempty(voxsiz);voxsiz=ones(1,3);end

assert(size(S,4)==32,'Visualization tools only implemented for 32 channel coils');

%VIRTUAL COIL PROFILES IN CENTRAL READOUT
Sr=S(ceil(end/2),:,:,:);
Src=compressCoils(Sr,size(S,4));%We assume there are 32 channels
Src=Src(:,:,:,1:16);
Src=squeeze(Src);
Src=permute(Src,[2 1 3]);
N=size(Src);
Src=reshape(Src,[N(1) N(2) 8 2]);
Src=permute(Src,[1 4 2 3]);
Src=reshape(Src,[N(1)*2 N(2)*8]);
FontSizeB=18;
FontSizeC=28;
figure
imshow(cat(1,abs(Src)/max(abs(Src(:))),(angle(Src)+pi)/(2*pi)),[])
if visSenPar~=2;title('Virtual sensitivity profiles in decreasing energy order','Interpreter','latex','FontSize',FontSizeC);end
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
if visSenPar==2
    t=1;
    path='/home/lcg13/Articulos/02EnProceso/Arxiv19-DISORDER/Figs/Fig04';
    export_fig(sprintf('%s/fig%02d.png',path,t));t=t+1;
    pause
end

%NORM OF 3D ROI VIRTUAL COIL PROFILES 
[Sc,~,~,D]=compressCoils(S,size(S,4));
D=sqrt(D);D=D/D(1);
D2=sqrt(multDimSum(abs(Sc).^2,1:3));D2=D2(:)';D2=D2/D2(1);
Dinf=multDimMax(abs(Sc),1:3);Dinf=Dinf(:)';Dinf=Dinf/Dinf(1);
figure
plot(D)
hold on
plot(D2)
plot(Dinf)
set(gca,'FontSize',FontSizeB)
xlabel('$c$','Interpreter','latex','FontSize',FontSizeC)
legend('$d^{1/2}$','$n^{[2]}_{\mathbf{U}^H\tilde{\bar{\mathbf{S}}}}$','$n^{[\infty]}_{a\mathbf{U}^H\tilde{\bar{\mathbf{S}}}}$','Interpreter','latex','FontSize',FontSizeC)
title('Normalized SVs of the decomposition ($d$) and norms of the components ($n$)','Interpreter','latex','FontSize',FontSizeC);
%axis([1 NP-1 -K(1)+1 K(1)-1])        
grid on
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
%export_fig(sprintf('%s/fig%02d.png',path,t));t=t+1;%In case we'd like to
%export the results

%SPECTRAL CONTENT OF 3D ROI VIRTUAL COIL PROFILES
[f,r]=radialPowerSpectralDensity(Sc,voxsiz);
str=cell(1,32);
figure
for s=1:32
    loglog(r,f(:,1,1,s))
    hold on
    str{s}=sprintf('$c=%s$',num2str(s));
end
grid on
set(gca,'FontSize',FontSizeB)
legend(str,'Interpreter','latex','FontSize',FontSizeB);
set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
ylabel('$|f|^2$','Interpreter','Latex','FontSize',FontSizeC);
xlabel('$|k|$ (mm$^{-1}$)','Interpreter','Latex','FontSize',FontSizeC);
title('Power spectral density','Interpreter','Latex','FontSize',FontSizeC);
pause

