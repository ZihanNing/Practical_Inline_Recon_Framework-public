function rec = Inline_AlignSENSE_recon_acs(twix)
%% 1. Generate the input of reconstruction (REC based on the twix-like structure
%%% RECOGNIZE SEQ TYPE
[seq_type,~] = recogSeqType(twix.hdr.sequenceParameters.sequence_type,...
twix.hdr.Meas.tScanningSequence,1);

%%% CONVERT TO REC STRUCTURE
fprintf('=========== Converting data to a reconstruction structure.\n');
fprintf('----------- Original raw read-in as converted twix-like data from ISMRMRD. \n');
if isequal(seq_type,'mege') 
    recon_flexFOV = [256 256 208];
    rec = twixlike2rec(twix,recon_flexFOV);
else
    rec = twixlike2rec(twix);
end
fprintf('Converted from twix-like data to input REC for following recon!\n');
fprintf('Acquisition name: %s.\n', rec.Names.Name);

% relief bufferData to save the memory 
% still keep the header info (twix.image,
% twix.bucket_header,twix.hdr,twix.sampling_description) to
% generate the header of the reconstructed image (see more details
% in referenceFromReconData_general.m
twix.data = []; 
if isfield(twix,'reference'); twix.reference = []; end
if isfield(twix,'noise'); twix.noise = []; end


%% 2. Estimate coil sensitivity map (load the exterREF)
% ZN: estimation method
solve_espirit = 1;
isFailed = 0;

fprintf('=========== Estimating coil sensitivity map using ACS line.\n'); 
sens_name = ['REF_DATA_ACS_',rec.Names.Name,'.mat']; % ZN: the ref file will contain seq name in the file name
nameRef = fullfile(rec.Names.pathOu, sens_name);

recS = rec; % ZN: borrow the rec structure
recS.y = rec.ACS;
recS.Enc.AcqVoxelSize = rec.Enc.UnderSampling.ACSVoxelSize;
recS.Par.Mine.APhiRec = rec.Par.Mine.APhiACS;
recS.Plan.Suff=''; recS.Plan.SuffOu='';
recS=solveSensit7T(recS); 
save(nameRef, 'recS','-v7.3');

%% 4. CALL RECONSTRUCTION
% pre-set parameters
if isfield(twix.hdr.acceleration,'type');accel_type = twix.hdr.acceleration.type;end % ZN: currently support GRAPPA (uniformed undersampling) & CAIPI undersampling
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
        NCon = size(rec.y,6);
        NAve = size(rec.y,5); 
        if isequal(seq_type,'mege') % For steady-state seq, we need to set the group motion states
            rec.Alg.parXT.sampleToGroup = round( rec.Par.Labels.TFEfactor/prod(rec.Enc.DISORDERInfo.tileSize));
        end
        recTemp = rec;
        % suppose single average for now
        p = []; % ZN: the image before moco
        x = []; % ZN: the image after moco
        for con = 1:NCon
            recTemp.Plan.Suff='_MotCorr';
            recTemp.Plan.SuffOu=sprintf('_Con%d',con);
            recTemp.y = dynInd(rec.y,con,6);
            if con==1 || seperate_recon
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
