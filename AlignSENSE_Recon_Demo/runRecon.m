
function [] = runRecon(idFile, idRef, idRefB, seq_type)

%% 0. Initialization
% List with all the names of the studies
studyNames;
B0In = fillCell(fileIn, '');
B1In = fillCell(fileIn, '');
isPT = fillCell(fileIn, 0);
supportReadout = fillCell(fileIn,[]);
resRec = fillCell(fileIn,[]);
if idRefB == 0
    refBIn{1}{1} = '';
    refBIn = fillCell(fileIn, '');
end
noiseFile = fillCell(fileIn,[]);

idPath= 1; 

supportReadoutRecon = fillCell(cell(1,length(idFile)),[]);
resRecRecon = fillCell(cell(1,length(idFile)),[]);

% Extract studies
[pathIn, fileIn, refIn, refBIn, B0In, B1In,...
fileUnique, refUnique, refBUnique, B0Unique, B1Unique, ...
isPTUnique,noiseFileUnique, supportReadoutUnique, resRecUnique, RDesiredUnique] ...
    =  extractStudies (pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath,idFile, idRef,[],[], isPT, noiseFile, supportReadout, resRec, RDesired);

% DATA CONVERSION PARAMETERS
estSensExplicit=0;
invDataExplicit = 0;

%% 1. Generate the input of reconstruction (REC based on the twix-like structure
%%% CONVERT TO REC STRUCTURE
for p=1:length(pathIn)  
    %BUILD THE ACQUISITION DATA - start with this if an acquisition file is also used as a reference, you don't read it ouy twice (dat2Se checks if rec structure exist whereas dat2Rec reads data out by default)
    for f=1:length(fileUnique{p})
        fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
        if (~exist(strcat(fileName,'.mat'),'file') || invDataExplicit) && ~strcmp(fileName,'')
            [rec, TW] = dat2Rec(fileName,supportReadoutUnique{p}{f},[],1,[],isPTUnique{p}{f},[],noiseFileUnique{p}{f},resRecUnique{p}{f},[],RDesiredUnique{p}{f});
        end                
    end
end
TW = [];

%% 2. Estimate coil sensitivity map 
% ZN: estimation method
solve_espirit = 1;
isFailed = 0;

fprintf('=========== Estimating coil sensitivity map using ACS line.\n'); 
sens_name = ['REF_DATA_ACS_',fileUnique{p}{f},'.mat']; % ZN: the ref file will contain seq name in the file name
nameRef = fullfile( rec.Names.pathOu, sens_name);

recS = rec; % ZN: borrow the rec structure
recS.y = rec.ACS;
recS.Enc.AcqVoxelSize = rec.Enc.UnderSampling.ACSVoxelSize;
recS.Par.Mine.APhiRec = rec.Par.Mine.APhiACS;
recS.Plan.Suff=''; recS.Plan.SuffOu='';
recS=solveSensit7T(recS); 
save(nameRef, 'recS','-v7.3');

%% 3. CALL RECONSTRUCTION
% pre-set parameters
accel_type = 'GRAPPA'; % ZN: currently support GRAPPA (uniformed undersampling) & CAIPI undersampling
rec.Alg.Gibbsring_flag=1; % ZN: do gibbs ringing reduction in matlab
rec.Alg.gibbsRinging = [0.4 0.4]; 
rec.Alg.CoilCompress_flag=0;  % ZN: do not use coil compression for MEGE sequence
% rec.Alg.CoilCompress_flag=config.CoilCompress_flag;  % ZN: do not use coil compression for MEGE sequence
seperate_recon = 1; %ZN: do the motion estimation for multiple contrasts seperately
% seperate_recon = config.seperate_recon; %ZN: do the motion estimation for multiple contrasts seperately
rec.Alg.disabledisplay = 1; % ZN: 1 - not display the image during the recon
rec.Alg.WriteSnapshots = 0; 

fprintf('=========== Reconstruction data.\n');
switch accel_type
    case 'GRAPPA'
        fprintf('The data is detected to be undersampled uniformly.\n');
        rec=assignSensitivities_bucket(rec, recS, solve_espirit);% Assign the sensitivities
        multiCon = size(rec.y,6)>1;
        multiEcho = size(rec.y,8)>1;
        if multiCon && multiEcho
            error('Cannot handling this currently.');
        else
            if multiCon
                NRep = size(rec.y,6);
            else
                NRep = size(rec.y,8);
            end
        end
        NAve = size(rec.y,5); 
        if isequal(seq_type,'mege') % For steady-state seq, we need to set the group motion states
            rec.Alg.parXT.sampleToGroup = round( rec.Par.Labels.TFEfactor/prod(rec.Enc.DISORDERInfo.tileSize));
        end
        recTemp = rec;
        % suppose single average for now
        p = []; % ZN: the image before moco
        x = []; % ZN: the image after moco
        for rep = 1:NRep
            recTemp.Plan.Suff='_MotCorr';
            if multiCon
                recTemp.Plan.SuffOu=sprintf('_Con%d',rep);
                recTemp.y = dynInd(rec.y,rep,6);
            else
                recTemp.Plan.SuffOu=sprintf('_Echo%d',rep);
                recTemp.y = dynInd(rec.y,rep,8);
            end
            if rep==1 || seperate_recon
               temp =  solveXT_gadgetron(recTemp);
               temp = gatherStruct(temp);
               p = cat(4,p,temp.x);
               rec.x = cat(4,x,temp.d);
               rec.T = temp.T;
               temp = [];
            else % Do not estimate the 
               temp = solveXT_gadgetron(recTemp,T);
               rec.p = cat(4,p,temp.x); % Aq, without MoCo
               rec.x = cat(4,x,temp.d); % Di, with MoCo
            end
        end
    otherwise % ZN: for CAIPI or other undersampling pattern
        fprintf('The data is detected to be undersampled nonuniformly.\n');
        NCon = size(rec.y,6);
        NAve = size(rec.y,5); 
        rec=assignSensitivities_bucket(rec, recS, solve_espirit,externalREF);% Assign the sensitivities
        flag_generateAq = 0; % ZN: for flexsampling pattern, the Aq image will need to be generated seperately
        [img_Di,img_Aq] = solveXT_gadgetron_FlexSamp(rec, flag_generateAq);
        if ~isempty(img_Aq);p = img_Aq; end
        rec.x=img_Di;
end

