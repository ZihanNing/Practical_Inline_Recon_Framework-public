addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'));

pathIn='/home/lcg13/Data/pnrawDe/ReconstructionsDebug07/2019_10_14/Ra_132330/Dy-Fu/';
fileIn{1}='ra_14102019_1755152_21_2_dhcp8sbfmrinorfnogradclearV4';
fileIn{2}='ra_14102019_1830541_25_2_dhcp8sbfmrinorfnogradnoch2preampclearV4';

compNoi=0;

for n=1:2%length(fileIn)
    for s=1:2
        load(strcat(pathIn,fileIn{n},sprintf('%03d',s),'.mat'));
        if s==1;recS{n}=rec;end
        if s==2;recS{n}.z=cat(5,recS{n}.z,rec.z);end
    end
    return
    recS{n}.z=permute(recS{n}.z,[1 2 8 5 4 3 6 7]);
    recS{n}.z=dynInd(recS{n}.z,recS{n}.Assign.z{8}+1,3);
    if compNoi
        xNoise=recS{n}.N(:,:); 
        xNoise=cat(1,real(xNoise),imag(xNoise));
        Ncov=cov(xNoise);
        Ncor=corrcoef(xNoise);    
        figure
        subtightplot(1,2,1)
        imshow(Ncov,[0 150])
        colormap(jet)
        colorbar
        title('Cov')
        subtightplot(1,2,2)
        imshow(Ncor)
        colormap(jet)    
        colorbar
        title('CorrCoef')
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
        export_fig(strcat(pathIn,'Snapshots',filesep,fileIn{n},'NoiseScan.png'))
        close all   
    
        xNoGrad=dynInd(recS{n}.z,size(recS{n}.z,4),4);
        NX=size(xNoGrad);
        xNoGrad=reshape(xNoGrad,[prod(NX(1:3)) prod(NX(4)) NX(5)]);    
        for l=1:size(xNoGrad,2)
            x=dynInd(xNoGrad,l,2);
            x=x(:,:);
            x=cat(1,real(x),imag(x));
            Ncov=cov(x);
            Ncor=corrcoef(x);    
            figure
            subtightplot(1,2,1)
            imshow(Ncov,[0 150])
            colormap(jet)
            colorbar
            title('Cov')
            subtightplot(1,2,2)
            imshow(Ncor)
            colormap(jet)    
            colorbar
            title('CorrCoef')
            set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
            export_fig(strcat(pathIn,'Snapshots',filesep,fileIn{n},'DataScanNoGrad.png'));
            close all
        end
    
        xNormal=dynInd(recS{n}.z,1:size(recS{n}.z,4)-1,4);
        NX=size(xNormal);
        xNormal=reshape(xNormal,[prod(NX(1:3)) prod(NX(4)) NX(5)]);
        for l=1:size(xNormal,2)
            x=dynInd(xNormal,l,2);
            x=x(:,:);
            x=cat(1,real(x),imag(x));
            Ncov=cov(x);
            Ncor=corrcoef(x);    
            figure
            subtightplot(1,2,1)
            imshow(Ncov,[0 150])
            colormap(jet)
            colorbar
            title('Cov')
            subtightplot(1,2,2)
            imshow(Ncor)
            colormap(jet)    
            colorbar
            title('CorrCoef')
            set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
            export_fig(strcat(pathIn,'Snapshots',filesep,fileIn{n},sprintf('DataScanStandard%03d.png',l)));
            close all
        end
    else
        xNormal=dynInd(recS{n}.z,1:size(recS{n}.z,4)/2,4);
        NX=size(xNormal);
        xNormal=reshape(xNormal,[prod(NX(1:4)) NX(5)]);
        NC=NX(5);
        figure
        for s=1:NC
            subtightplot(NC,1,s)
            plot(abs(xNormal(:,s)))
        end
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
        export_fig(strcat(pathIn,'Snapshots',filesep,fileIn{n},'ProfileDataScanStandard.png'));
    end
end

