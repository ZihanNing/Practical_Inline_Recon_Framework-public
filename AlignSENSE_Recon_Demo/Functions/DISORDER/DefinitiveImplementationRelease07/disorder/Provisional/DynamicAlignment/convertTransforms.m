%%path{1}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/17Nov17';
%%path{2}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/13Nov17';
%%path{3}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/10Nov17';
%%path{4}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/10Nov17b';
%%path{5}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/6Nov17b';
%%path{6}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/3Nov17';%%%
path{1}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/1Nov17';%%%

%%path{1}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/13Oct17b';
%%path{2}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/9Oct17b';
%%path{3}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/18Sept17b';
%%path{4}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/8Sept17b';
%%path{5}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/8Sept17';
%%path{6}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/6Sept17b';%%%
%%path{7}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/21Aug17';%%%
%%path{8}='Data/rawSource/dhcp-share-kcl/fetal_workup/FMRI/22Nov17-2';

for n=1:length(path)
    n
    pathF=strcat('/home/lcg13/',path{n},'/Dy-Fu2');
    cd(pathF);    
    fileTr=dir('*_Tr.mat');
    for m=1:length(fileTr)
        load(fileTr(m).name);
        T=gather(T);
        save(fileTr(m).name,'T');
    end
    cd('/home/lcg13/Work/DefinitiveImplementationDebug05');  
end
