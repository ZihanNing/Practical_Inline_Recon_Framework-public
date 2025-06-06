
function [pathIn, fileList,       refList,       refBList,      B0List,       B1List,...
                  fileListUnique, refListUnique, refBListUnique, B0ListUnique, B1ListUnique,...
                  isPTListUnique, noiseFileListUnique, supportReadoutListUnique, resRecListUnique, RDesiredListUnique] ...
                  = extractStudies(pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath, idFile, idRef, idB0, idB1, isPT, noiseFile, supportReadout, resRec, RDesired)
                           
%EXTRACTSTUDIES extractst the filenames for the differnt acquisitions, based on a selection provided by the used.
%   []=EXTRACTSTUDIES()
% 

assert(length(idFile)==length(idRef),'extractStudies:: idRef and idFile should have the same size.')
if nargin < 10 || isempty(idB0); idB0 = onesL(idRef);end
if nargin < 11 || isempty(idB1); idB1 = onesL(idRef);end
if nargin < 12 || isempty(isPT); isPT = fillCell(fileIn,0);end
if nargin < 13 || isempty(noiseFile); noiseFile = fillCell(fileIn,[]);end
if nargin < 14 || isempty(supportReadout); supportReadout = fillCell(fileIn,[]);end
if nargin < 15 || isempty(resRec); resRec = fillCell(fileIn,[]);end
if nargin < 16 || isempty(RDesired); RDesired = fillCell(fileIn,[]);end

%%% EXTRACT PATH
if ~isempty(idPath)
    idPath(idPath>length(pathIn))=[];
    pathIn=pathIn(idPath);
    fileIn=fileIn(idPath);
    refIn=refIn(idPath);
    refBIn=refBIn(idPath);
    B0In=B0In(idPath);
    B1In=B1In(idPath);
end

%%% Make list with all acquisitions to be processed in dat2Rec
fileListUnique=[];
refListUnique=[];
refBListUnique=[];
B0ListUnique=[];
B1ListUnique=[];
isPTListUnique=[];
noiseFileListUnique=[];
supportReadoutListUnique=[];
resRecListUnique=[];
RDesiredListUnique=[];

for p =1:length(pathIn)
    
    %%% MAKING LIST OF UNIQUE ACQUISITIONS TO PROCESS
    if iscell(idFile)
        fileListIdx = [];
        for i =1:length(idFile); fileListIdx = cat(2 ,fileListIdx, idFile{i});end
        fileListIdx = unique(fileListIdx);
    else 
        fileListIdx = unique(idFile);
    end
    for i=find(fileListIdx)
        fileListUnique{p}{i} = fileIn{p}{fileListIdx(i)};
        isPTListUnique{p}{i} = isPT{p}{fileListIdx(i)};%Only the unique one since 
        noiseFileListUnique{p}{i} = noiseFile{p}{fileListIdx(i)};
        supportReadoutListUnique{p}{i} = supportReadout{p}{fileListIdx(i)};
        resRecListUnique{p}{i} = resRec{p}{fileListIdx(i)};
        RDesiredListUnique{p}{i} = RDesired{p}{fileListIdx(i)};
    end%create cell array with names
 
    %%% ACQUISITIONS
    fileList=[];
    if ~iscell(idFile)
        idFile(idFile>length(fileIn{p}))=[];%if idFile is larger than actuall # of acquisitions, leave those out     
        fileList{p}=fileIn{p}(idFile);
    else
        for i=1:length(idFile)
            fileList{p}{i} = fileIn{p}(idFile{i});
        end
    end

    %%% REFERENCE
    if ~iscell(idRef)
        idRef(idRef>length(refIn{p}))=[];%if idFile is larger than actuall # of acquisitions, leave those out 
        refList{p}=refIn{p}(idRef);
        refBList{p}=refBIn{p}(idRef);
    else
        error('extractStudies:: idRef should not be cell array.')
    end
    refListIdx = unique(idRef);
    refBListIdx = unique(idRef);
    for i=find(refListIdx); refListUnique{p}{i} = refIn{p}{refListIdx(i)};refBListUnique{p}{i} = refBIn{p}{refBListIdx(i)};end%create cell array with names
     
    %%% B0 MAPS
    if ~iscell(idB0)
        idB0(idB0>length(B0In{p}))=[];%if idFile is larger than actuall # of acquisitions, leave those out   
        B0List{p}=B0In{p}(idB0);
    else
        error('extractStudies:: idB0 should not be cell array.')
    end
    B0ListIdx = unique(idB0);
    for i=find(B0ListIdx); B0ListUnique{p}{i} = B0In{p}{B0ListIdx(i)};end%create cell array with names
    
    %%% B1 MAPS
    if ~iscell(idB1)
        idB1(idB1>length(B1In{p}))=[];%if idFile is larger than actuall # of acquisitions, leave those out  WRONG, could be multiple recons on same data, so CHANGE   
        B1List{p}=B1In{p}(idB1);
    else
        error('extractStudies:: idB1 should not be cell array.')
    end
    B1ListIdx = unique(idB1);
    for i=find(B1ListIdx); B1ListUnique{p}{i} = B1In{p}{B1ListIdx(i)};end%create cell array with names
end