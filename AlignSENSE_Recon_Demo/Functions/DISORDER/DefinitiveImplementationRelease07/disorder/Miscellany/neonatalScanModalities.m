
addpath(genpath(fileparts(mfilename('fullpath'))));
[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
caseFetal=0;
[paths,pathid]=studiesDetection(1+caseFetal);
%paths=paths(1:2);

modalities={'pilotV4','refneoheadV4','b0mapV4','b1mapclearV4','t2tsesenseV4','mpragesenseV4','sefmriclearV4','sbfmriclearV4','mbfmriV4','sbfmrirepclearV4','sbdtisenseV4','mbdti2senseV4','b0mapshimV4',...
    't1tseirsenseV4',... %Standard scans up to this point
    'mbdtisenseV4','b1mapdefaultclearV4','highrest1volclearV4','sbfmrite451700clearV4','mb9fte451700experiment1V4','sbfmrioldV4','mbfmrioldV4','dynaneonate1senseV4','dti32senseV4',...
    'mprage2clearV4','sbfmri2clearV4','mbdti3senseV4','arteriogramsenseV4','pca06x4ga40clearV4','pca06x3ga40se2senseV4','b0mappbV4','sefmrirepclearV4','sbdtirepsenseV4','sefmrirclearV4','sbclearV4',...
    'mb3grte98clearV4','mb2grte115clearV4','sefmripaclearV4','sefmriapclearV4','sefmriappclearV4','dw24senseV4','dw6natlrsenseV4','dw6natrlsenseV4','dw6natpasenseV4','dw6natapsenseV4','dw12natrl12senseV4',...
    'dw12natap12senseV4','eyesbsesenseV4','eyesbgeclearV4','press55dynamic42pbvolumefirstordershimmingV4','press144dynamic135pbvolumefirstordershimmingV4','press288dynamic270pbvolumefirstordershimmingV4',...    
    'mbdti3000senseV4','sbdtiisosenseV4','mbdti3000isosenseV4','angiogramsenseV4','dsbap1000defprosenseV4','dsbap1000optteb6senseV4','mb4ap1000senseV4','mbdtiisosenseV4','dsbap1000restsenseV4',...
    'mb4ap1000restsenseV4','survey32chheadcoilV4','swite25senseV4','swite25r050senseV4','svpress144bgV4','svpress288bgV4','dwisenseV4','neoquickadcsenseV4','pilotqbcV4','pilotquietV4','pilotpredticlearV4',...
    'neckplanscansensesenseV4','sbhighresclearV4','t2w3ddrive32chshcsenseV4','pcaV4','angiocowsenseV4',...%Using the dHCP patch up to this point
    'wipven3dpcasenseV4','wiparteriogramsenseV4','wipswite25r050senseV4','wipsensehiresfusenseV4','wipneckplanscansenseV4','ven3dpcasenseV4','sensehiresfusenseV4','t2senseV4'};
typ=cell(1,length(modalities));
for n=1:length(modalities);typ{n}='uint8';end
variableNames=['scanNo','subId','sesId',modalities];
variableTypes=['uint16','string','string',typ];
M=table('Size', [length(paths) length(variableNames)],'VariableNames',variableNames,'VariableTypes',variableTypes);

N=cell(length(paths),length(variableNames));

for s=1:length(paths)
    fprintf('Study %d of %d\n',s,length(paths));
    pathIn=rawFolderDetection(paths{s});
    dirContents=dir(fullfile(pathIn,'*.lab'));
    dirContents={dirContents.name};
    a=zeros(1,length(modalities));
    for l=1:length(dirContents)
        foundMod=0;
        for m=1:length(modalities)
            if contains(dirContents{l},modalities{m})
                a(m)=a(m)+1;
                foundMod=1;
                break
            end
        end   
        if ~foundMod
            fprintf('Unknown modality:\n');
            fprintf('%s\n',pathIn);
            fprintf('%s\n',dirContents{l});
            return
        end
    end
    M.scanNo(s)=s;
    for m=1:2;M.(variableNames{m+1})(s)=pathid{s}{m};end
    for m=1:length(modalities);M.(modalities{m})(s)=a(m);end
end

outPath='/home/lcg13/Data/pnrawDe/ReconstructionsRelease03/INF/';
fileNameOu=sprintf('%s/scanAttemptedModalities.csv',outPath);
writetable(M,fileNameOu);
