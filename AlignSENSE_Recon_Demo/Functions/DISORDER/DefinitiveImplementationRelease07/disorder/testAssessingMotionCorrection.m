
addpath(genpath(fileparts(mfilename('fullpath'))));


inpPath='/home/lcg13/Data/pnrawDe';
outFol='TS';

load(sprintf('%s/%s/anonKeys.mat',inpPath,outFol))
sTInfo=[];
sedStr={'Sedated','Not sedated'};
motStr={'Motion corrected','Not motion corrected'};
for n=1:size(anonKeys,1)
    sTInfoC{1}=sprintf('%04d',n);
    sTInfoC{2}=sedStr{floor((anonKeys(n,1)+1)/2)};
    sTInfoC{3}=motStr{mod(anonKeys(n,1),2)+1};
    sTInfoC{4}=pathr{floor((anonKeys(n,1)+1)/2)}{anonKeys(n,2)}{1};
    sTInfo=cat(1,sTInfo,sTInfoC);
end
    
%sTInfo=sTInfo(:,[2 3 6 7 5 4 1]);
sTInfoTable=table(sTInfo);
writetable(sTInfoTable,'tableSedationVsMotionCorrection.xls');
    
    
%    anonKeys(n,1)
    
%end


return



[user,versCode,versSave,pathRF,pathSt,pathPr]=versRecCode;
motionCorrectionAssess;

for n=1:length(path)
    cr=1;    
    for c=1:length(path{n})        
        if length(path{n}{c})==2
            dData=sprintf('%s/ReconstructionsRelease03/sub-%s/ses-%s/T2MS',inpPath,path{n}{c}{1},path{n}{c}{2});
            structFile=dir(sprintf('%s/*Co3DOutSVRMot.nii',dData));           
            structFile={structFile.name};           
            if ~isempty(structFile);dFile=structFile{1};else dFile=[];end
        elseif n==1
            dData=sprintf('%s/ChrisRelease03/%s/T2MS',inpPath,path{n}{c});
            dFile='reconT2MSCo.nii';
            if ~exist(sprintf('%s/%s',dData,dFile),'file');dFile=[];end          
        else
            dData=sprintf('%s/ReconstructionsRelease03/sub-%s',inpPath,path{n}{c});
            if exist(dData,'dir')
                structPath = dir(dData);%Returns all the files and folders in the directory
                structPath(ismember( {structPath.name}, {'.', '..'}))=[];
                structPath=structPath([structPath.isdir]);                
                if ~isempty(structPath);dData=sprintf('%s/%s/T2MS',dData,structPath(1).name);else dData=[];end
            else dData=[];
            end
            if ~isempty(dData)
                structFile=dir(sprintf('%s/*Co3DOutSVRMot.nii',dData));           
                structFile={structFile.name};           
                if ~isempty(structFile);dFile=structFile{1};else dFile=[];end
            else dFile=[];
            end
        end
        if ~isempty(dFile)
            if length(path{n}{c})==2 || n==2;structFile=dir(sprintf('%s/*_UnMot.nii',dData));else structFile=dir(sprintf('%s/*_Un.nii',dData));end
            structFile={structFile.name};
            if length(structFile)==2
                pathr{n}{cr}{1}=dData;
                pathr{n}{cr}{2}=dFile;
                for s=1:2;pathr{n}{cr}{2+s}=structFile{s};end
            end
            cr=cr+1;
        end
    end
end

NCasesTotal=(length(pathr{1})+length(pathr{2}))*3;
NCasPerm=randperm(NCasesTotal);

cr=1;
outFol='TestsSedation';
anonKeys=zeros(NCasesTotal,2);%1->Sedated non motion corrected/2->Sedated motion corrected/3->Non sedated non motion corrected, 4->Non sedated motion corrected / Second element refers to case in the 
for n=1:length(pathr)
    for c=1:length(pathr{n})
        for s=1:3
            labId=NCasPerm(cr);
            outFile=sprintf('%s/%s/%04d.nii',inpPath,outFol,labId);
            inpFile=sprintf('%s/%s',pathr{n}{c}{1},pathr{n}{c}{1+s});
            copyfile(inpFile,outFile);
            if s==1 && n==1
                anonKeys(labId,1)=2;%Sedated motion corrected
            elseif s==1 && n==2
                anonKeys(labId,1)=4;%Non sedated motion corrected
            elseif n==1
                anonKeys(labId,1)=1;%Sedated non motion corrected
            else
                anonKeys(labId,1)=3;%Non sedated non motion corrected
            end
            anonKeys(labId,2)=c;
            cr=cr+1;
        end
    end
end

save(sprintf('%s/%s/anonKeys.mat',inpPath,outFol),'anonKeys','pathr');

%load(sprintf('%s/%s/anonKeys.mat',inpPath,outFol))

return





inpRootPath=[];%'/home/lcg13/Data/18ExperimentsTracking';%[];%inpRootPath='/home/lcg13/Data/rawSource';%Modify to specify the input (raw) path
shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Neonatal';
outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/Reference';%/home/lcg13/Data/18ExperimentsTracking';%[];%'/home/lcg13/Data/rawDestin/ReconstructionsDebug06';%outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug05';%Modify to specify the output (reconstructions) path
if caseFetal
    shaRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug07/Data/Fetal';
    outRootPath='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/ReferenceFetal';
end
modality=[2:4 9];%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
series=[];%[];%[];%[];%3;%1:100;%2:10;%3:10;%[];%[];%All series. Modify to specify a particular series (with indexes as row position in prot.txt)

NCori=1;
NCend=length(paths);

%for n=NCori:NCend
%    fprintf('\n\nCase %d of %d: %s\n',n-NCori+1,NCend,paths{n});
%    parseLabFolder(paths{n},modality,outRootPath,inpRootPath);%PARSING RAW REQUIRES RECONFRAME LICENSE
%end

%return

for n=NCori:NCend%468%[1 425 NCend]%:1%length(paths)
    modality=2:4;%2:4;%[];%[];%[];%All modalities. Modify to specify a particular modality (see list in Control/modalList.m)
    tsta=tic;
    fprintf('\n\nCase %d of %d: %s',n-NCori+1,NCend,paths{n});
    fRec=protoPipeline(paths{n},modality,outRootPath,inpRootPath,series);

    modality=[];  
    fileNameConversion(paths{n},pathid{n},modality,outRootPath,inpRootPath,shaRootPath,series);
    tend=toc(tsta);fprintf('Time: %.3f s\n',tend);
end