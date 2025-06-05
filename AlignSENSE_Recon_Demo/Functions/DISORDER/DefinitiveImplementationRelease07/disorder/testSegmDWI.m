
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;

paths=studiesDetection(2);

suffIn{1}='AveBMax';
suffIn{2}='AveBMaxMaskFirst';
gpu=1;

%SEGMENTATION ISSUES:
%Case 26 leaks a bit
%Case 46 slightly undermasked
%Case 48 leaks a bit
%Case 85 really challenging

%Check for rings:2018_07_04/SA_24430(82)
%paths{1}='2018_07_04/SA_24430';

for p=114:length(paths)
    pathRIn=strcat('/home/lcg13/Data/pnrawDe/ReconstructionsRelease06SVR/',paths{p},'/Dy-Di/');
    pathROu=strcat('/home/lcg13/Data/pnrawDe/ReconstructionsRelease06SVR/',paths{p},'/Dy-Di/');
    structPath = dir(pathRIn);%Returns all the files and folders in the directory
    structPath(ismember( {structPath.name}, {'.', '..'}))=[];
    structFil=structPath(~[structPath.isdir]);

    c=0;
    for n=1:length(structFil)        
        if length(structFil(n).name)>12 && strcmp(structFil(n).name(end-11:end),'_AveBMax.nii')
            c=1;
            candF=structFil(n).name;           
            break;
        end
    end
    if c>0        
        fileRECONIn=strcat(pathRIn,candF(1:end-12));
        fileRECONOu=strcat(pathROu,candF(1:end-12));
        fprintf(sprintf('%s\n',fileRECONIn));
        tic
        [y,MS,MT]=readNII(fileRECONIn,suffIn,gpu);
        toc

        x=y{1};
        Mla=y{2};

        %ELLIPSOID FIT
        tic
        [par,Mlb]=ellipsoidFromImage(Mla,[],0);
        toc
        
        tic
        K=32;
        S=128;
        %rRange=[5 30];
        rRange=[5 30];%RECENT CHANGE
        like=[0 -1 -1 0];
        %pri=[0 0 0 0 5 2];%RECENT CHANGE
        pri=[0 0 0 0 20 100];%RECENT CHANGE
        grTh=[0.55 0.95;%Threshold (relative to the  maximum gradient) for gradient feature extraction
             0.05 0.95];%Percentile for gradient intensity feature extraction
        extK=[2 2];
        diUse=2;
        kappa=[0.1 0.1];
        eroDilaFactor=[3 4];
        [Mlc,par,F,xF]=brainSegmentation(x,par(1:3),rRange,K,S,mean(par(4:6)),like,pri,diUse,kappa,extK,grTh,eroDilaFactor);%diUse,kappa,extK,grTh
        z{1}=Mlc;
        %z{2}=F{1};
        %z{3}=F{4};
        %z{4}=xF{1};
        suff{1}='AveBMaxMaskThird';
        %suff{2}='Int';
        %suff{3}='Gra';
        %suff{4}='Filt';
        toc
        
        tic
        NS=length(suff);
        MSV=cell(1,NS);
        MTV=cell(1,NS);
        for n=1:NS;MSV{n}=MS{1};MTV{n}=MT{1};end
        writeNII(fileRECONOu,suff,z,MSV,MTV);     
        
        %return
    end
end
return


[ME,sphcen0,sphrad]=ellipsoidalHoughTransform(feat,voxsiz,[20 60]);
rec.Par.Mine.sphcen0=sphcen0;
rec.Par.Mine.sphrad=sphrad;
if N(4)>10
    writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'ME','feat'},{ME,feat});
    %1
    %pause
end




