function Par=lab2Mat(Param,modal,modals,curStud)

%LAB2MAT   Converts a ReconFrame (RF) metadata object into a Matlab 
%structure
%   PAR=LAB2MAT(PARAM,MODAL,{MODALS},{CURSTUD})
%   * PARAM is a ReconFrame object
%   * MODAL is a modality to be assigned to the metadata
%   * {MODALS} is the set of potential modalities available
%   * {CURSTUD} is the information for specific studies
%   ** PAR is the returned matlab structure
%

%SET DEFAULT VALUES AND CHECK FOR ERRORS
if nargin<3;modals=[];end
if nargin<4 || isempty(curStud);curStud.Stud='None';end

if isempty(modals) || ~isnumeric(modals);modals=modal2Numbe(modals);end
assert(numel(modal)==1,'Only a single modality can be assigned');
if ~isnumeric(modal);modal=modal2Numbe(modal);end

%INITIALIZE THE OUTPUT STRUCTURE AND CHECK FOR MODALITY IN SET OF POTENTIAL
%MODALITIES, ALSO WE INCLUDE HERE THE ADHOCARRAY
Par=[];
if ismember(modal,modals);Par.Mine.Modal=modal;else return;end
if isfield(Param.Labels,'AdHocArray');Par.Mine.AdHocArray=Param.Labels.AdHocArray;
elseif Param.IsParameter('MP_RECFRAME_ad_hoc');Par.Mine.AdHocArray=Param.GetValue('MP_RECFRAME_ad_hoc');
else Par.Mine.AdHocArray=[];
end

[path,name]=fileparts(Param.Filename.Data);
path=strsplit(path,filesep);
path=strcat(path{end-1},filesep,path{end});
if strcmp(path,'2019_07_08/O__111430') && strcmp(name,'oMu_08072019_1823268_5_2_mv3dmrisb2mmsenseV4');Par.Mine.AdHocArray(1:128)=0;Par.Mine.AdHocArray(2)=20;
elseif strcmp(path,'2019_07_08/O__111430') && strcmp(name,'oMu_08072019_1827051_6_2_mv3dmri2mmsenseV4');Par.Mine.AdHocArray(1:128)=0;Par.Mine.AdHocArray(1:8)=[101 50 3 3 44 0 2 1];
elseif strcmp(path,'2019_11_01/Ja_136730') && strcmp(name,'ja_01112019_1724132_22_2_fiss3mb3btfebhfissprofclearmo2V4');Par.Mine.AdHocArray(4)=5;
elseif strcmp(path,'2019_11_01/Ja_136730') && strcmp(name,'ja_01112019_1738273_26_2_fiss3mb3btfebhclearmo2V4');Par.Mine.AdHocArray(4)=5;
end 

%BRUTE-FORCE CONVERSION OF RF OBJECT INTO MATLAB STRUCTURE
if isobject(Param)
    lev1Fields=fieldnames(Param);
    for n=1:length(lev1Fields)
        curFieldLev1=Param.(lev1Fields{n});
        if isobject(curFieldLev1)
            lev2Fields=fieldnames(curFieldLev1);
            for m=1:length(lev2Fields)                
                Par.(lev1Fields{n}).(lev2Fields{m})=Param.(lev1Fields{n}).(lev2Fields{m});
            end
        else
            Par.(lev1Fields{n})=Param.(lev1Fields{n});
        end
    end
else
    error('The intended RF object is actually not an object');
end

%READ PES, B-VALS AND ORIENTATIONS FROM FILE IF NECESSARY
%THE FOLLOWING FILES ARE / WILL BE CONTEMPLATED:
%AP_PA.txt: directions of auxiliary spin echo data used for fMRI distortion correction in neonatal dHCP, 8 vols, not sure why an additional line is present
%dhcpneoSB.txt: directions of SB DWI data in neonatal dHCP, 8 vols
%dhcp300_f.txt: directions of DWI data in neonatal dHCP, 300 vols
%dhcp_fetal142_c.txt: directions of DWI data in current fetal dHCP with two b0s at the beginning, 142 vols
%dhcp_fetal140.txt: directions of DWI data in current fetal dHCP, 140 vols
%fet_bval_1000i.txt: directions of DWI data in reduced fetal dHPC, 29 vols
%dhcp_fetal90.txt: directions of DWI data in deliberately slow fetal dHCP, 90 vols
%dw_scheme100.txt: directions of DWI in infant studies, 100 vols
%dw_scheme150.txt: directions of DWI in infant studies, 150 vols
%dw_scheme_b10000.txt: directions of DWI in adult tests, 550 vols
%restart_dirs_15_1150210_1134_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 193 vols
%restart_dirs_16_1150210_1152_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 113 vols
%restart_dirs_16_co_07052015_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 157 vols
%restart_dirs_16_na_01042015_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 156 vols
%restart_dirs_17_co_07052015_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 149 vols
%restart_dirs_17_na_01042015_dhcp300_f.txt: directions of DWI for corresponding neonatal dHCP study, 156 vols
Par.Mine.diInfo=[];Par.Mine.curStud=curStud;
if ~strcmp(curStud.Stud,'None') && Par.Mine.Modal==10
    if strcmp(curStud.Stud,'dHCPTodBra') && Par.Mine.AdHocArray(1)==0 && strcmp(curStud.IssuId,'MismatchedReference');direFile='toddlerSinglePE.txt';%NOT SURE THIS WAS WORKING WELL
    elseif (((strcmp(curStud.Stud,'dHCPNeoPha') || strcmp(curStud.Stud,'dHCPNeoBra') || strcmp(curStud.Stud,'dHCPTodPha') || strcmp(curStud.Stud,'dHCPTodBra')) && Par.Labels.PartialFourierFactors(1,2)~=1) || strcmp(curStud.Stud,'dHCPAduBra') || strcmp(curStud.Stud,'dHCPEpiPha') || strcmp(curStud.Stud,'dHCPEpiBra')) && Par.Mine.AdHocArray(1)==0;direFile='dhcpneoSB.txt';
    elseif (strcmp(curStud.Stud,'dHCPNeoPha') || strcmp(curStud.Stud,'dHCPNeoBra') || strcmp(curStud.Stud,'dHCPTodPha') || strcmp(curStud.Stud,'dHCPTodBra')) && Par.Labels.PartialFourierFactors(1,2)==1;direFile='AP_PA.txt';
    elseif (strcmp(curStud.Stud,'dHCPNeoPha') || strcmp(curStud.Stud,'dHCPNeoBra')) && Par.Mine.AdHocArray(1)==101;direFile='dhcp300_f.txt';
    elseif (strcmp(curStud.Stud,'dHCPTodPha') || strcmp(curStud.Stud,'dHCPTodBra')) && strcmp(curStud.IssuId,'None') && Par.Mine.AdHocArray(1)==101;direFile='anx_noddi_ap.txt';%direFile='anx_noddi.txt';%direFile='dw_anxiety.txt';
    elseif strcmp(curStud.Stud,'dHCPAduBra') && Par.Mine.AdHocArray(1)==101 && (strcmp(curStud.IssuId,'None') || strcmp(curStud.IssuId,'DiscardDWI'));direFile='dw_scheme_b10000Corrected.txt';
    elseif strcmp(curStud.Stud,'dHCPAduBra') && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'None1');direFile='scheme450short.txt';
    elseif strcmp(curStud.Stud,'dHCPAduBra') && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'None2');direFile='scheme450.txt';
    elseif strcmp(curStud.Stud,'dHCPTodBra') && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'MismatchedReference');direFile='dw_scheme150SinglePE.txt';%THIS HAS BEEN USED TO TEST WHAT WAS GOING ON WITH CASE 2018_06_23/WI_22230--IT IS REALLY NOT MATCHING DW_SCHEME150!!
    elseif strcmp(curStud.Stud,'dHCPTodBra') && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'OldData');direFile='dw_scheme150.txt';
    elseif strcmp(curStud.Stud,'dHCPTodBra') && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'100Dir');direFile='dw_scheme100.txt';
    %elseif (strcmp(curStud.Stud,'dHCPEpiBra') || strcmp(curStud.Stud,'dHCPEpiPha')) && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'None');direFile='dw_scheme.txt';
    elseif (strcmp(curStud.Stud,'dHCPEpiBra') || strcmp(curStud.Stud,'dHCPEpiPha')) && Par.Mine.AdHocArray(1)==101 && strcmp(curStud.IssuId,'None');direFile='rice.txt';%Epilepsy
    else fprintf('Undefined study type %s',curStud.Stud);
    end
    [~,~,~,~,pathSt]=versRecCode;
    diFile=fullfile(pathSt,'DynamicFiles',direFile);
    if exist(diFile,'file');Par.Mine.diInfo=load(diFile);end
end

%CHECK CONSISTENCY AND PERFORM HEADER FEEDING
assert(isfield(Par,'Labels'),'No Labels for file %s',Par.Filename.Parameter);
assert(isfield(Par.Labels,'Index'),'No Indexes for file %s',Par.Filename.Parameter);
%feedHeader(Par,1);%TO CHECK CONSISTENCY IN CASE PROBLEMS ARE ENCOUNTERED
if modal==10;Par=feedHeader(Par);end%WE ONLY FEED IT FOR DIFFUSION AT THE MOMENT

%REDUCE THE PROFILE INFORMATION BY COMPRESSING THE LABELS AMONG THE DIFFERENT COILS
%CHECK THAT THAT THE COMPRESSION IS POSSIBLE---THIS IS NOW DISABLED
Index=Par.Labels.Index;
Par2Read=Par.Parameter2Read;
NC=length(Par.Parameter2Read.chan);%Number of surface coils
Par.Mine.SurfaceCoils=NC;
if length(Par2Read.loca)==2 && length(Par2Read.kz)>1;NC=find(diff(single(Index.loca))==-1,1);end%Most likely a reference scan->Number of surface + body coils
if NC==0;Par=[];return;end
Par.Mine.BodyCoils=NC-Par.Mine.SurfaceCoils;
% if Par.Mine.Modal~=2 && Par.Mine.Modal~=3 && Par.Mine.Modal~=5 && Par.Mine.Modal~=6 && Par.Mine.Modal~=7
%     NP=length(Index.chan);
%     indexNames=fieldnames(Index);
%     for m=1:length(indexNames)
%         subIndexV=Index.(indexNames{m});
%         assert(length(subIndexV)==NP,'Size of %s (%d) does not match size of chan (%d)',indexNames{m},length(subIndexV),NP);
%         assert(mod(NP,NC)==0,'Number of profiles (%d) divided by number of coils (%d) is not an integer for the %s in file %s',NP,NC,indexNames{m},Par.Filename.Parameter);
%         subIndexM=reshape(subIndexV,[NC NP/NC]);
%         if strcmp(indexNames{m},'chan');subIndexM=subIndexM';end
%         %This is for coil data, not used as we are not compressing the
%         %information. However it is left to detect potential inconsistencies in
%         %the future. It was OK for MRecon-3.0.460 but in MRecon-3.0.515 some
%         %new profiles are added after the noise acquisition. They are
%         %identified as typ 4, correspond to loca 0, body coil chan ids and 
%         %chan_grp ids
%         %if ~strcmp(indexNames{m},'offset') && (Par.Mine.Modal~=2 || (~strcmp(indexNames{m},'loca') && ~strcmp(indexNames{m},'rtop') && ~strcmp(indexNames{m},'random_phase') && ~strcmp(indexNames{m},'pda_fac') && ~strcmp(indexNames{m},'chan_grp'))) 
%         if ~strcmp(indexNames{m},'offset') && ~strcmp(indexNames{m},'pda_fac')%This line would replace the
%      NC   %previous one if these lines are not run for coil data
%             assert(size(unique(subIndexM,'rows'),1)==1,'The %s are not ordered independently from the channels for file %s',indexNames{m},Par.Filename.Parameter);
%         end
%     end
% end

% %PERFORM THE COMPRESSION OF INDEXES. Note the offsets are not compressed,
% %perhaps because they relate with the lenght of the readout and
% %this is bigger for noise samples. Considering the code for checking
% %compression introduced before, despite the coil metadata has not been
% %compressed, it could be, at least for those indexes other than the loca, 
% %rtop, random_phase, pda_fact and chan_grp. Thus, there could be room for
% %compression on both sides, but probably only providing marginal benefits.
% %Further tests have shown that compression is not possible either for
% %structural data, where we have found that there is repeated acquisition of
% %typ 3 data over what seems to be the body coils. 3000 samples where
% %acquired for each of the body coils for a total sampling space of 72
% %k-space lines x 125 slices; they correspond to the repeated acquisition of
% %the central line (1/3x the number of PE lines) for each slice. As for the 
% %T1MS these number are 1750 / 70 / 125, i.e., 1/5x the number of PE lines).
% %The T2MS has 6 shots while the T1MS has 10 shots, so it seems they have
% %been acquired for half the number of shots. Similar behaviour for 
% %MPRAGE, but only 2 profiles of typ 3 acquired (one for each channel of 
% %the body coil probably?). Similar for B0 data in '2017_01_18/ge_193200',
% %with the number of typ 3 profiles being a multiple of the number of
% %slices. For SE the pda_fac are not ordered independently from the 
% %channels, so they need to be stored as the offsets. As for the T1MS
% if Par.Mine.Modal~=2 && Par.Mine.Modal~=3 && Par.Mine.Modal~=5 && Par.Mine.Modal~=6 && Par.Mine.Modal~=7
%     for m=1:length(indexNames)
%         subIndexV=Par.Labels.Index.(indexNames{m});
%         if ~strcmp(indexNames{m},'offset') && ~strcmp(indexNames{m},'pda_fac')            
%             if strcmp(indexNames{m},'chan')                
%                 Par.Labels.Index.(indexNames{m})=subIndexV(1:NC);
%             else
%                 Par.Labels.Index.(indexNames{m})=subIndexV(1:NC:end);
%             end
%         end
%     end
% end

% %DUMPING THE STRUCTURE
warning('off','MATLAB:structOnObject');
ParStruct=struct(Param);
Par.GoalC=[];
if isfield(ParStruct,'Par40') && ~isempty(ParStruct.Par40);
    %WRITING GOAL-C PARAMETERS OF MOST INTEREST (I HAVEN'T BEEN ABLE TO FIND THEM LATER...)
    %DUAL SPIN-GRADIENT ECHO
    if size(Param.Encoding.YRange,1)==2 && Param.Encoding.NrMixes==1      
        Par.Mine.Spin2Echo=Param.IsObject('RF`ME');%Indicates whether the second echo is a spin echo
        if ~Par.Mine.Spin2Echo && Param.IsObject('GR`pyr_ME');Par.Mine.StrFactorMax=Param.GetValue('GR`pyr_ME:str_factor_max');end
    end
    %SHIM BOX CENTER
    if Param.IsParameter('EX_AS_vols_ap_offcentres') && Param.IsParameter('UGN1_AS_vols') && Param.GetValue('UGN1_AS_vols')==1;Par.Mine.ShimBoxCenter=[Param.GetValue('EX_AS_vols_lr_offcentres') Param.GetValue('EX_AS_vols_ap_offcentres') Param.GetValue('EX_AS_vols_fh_offcentres')];end    
    %ECHO SPACING
    if Param.IsParameter('UGN3_ACQ_epi_es');Par.Mine.ES=Param.GetValue('UGN3_ACQ_epi_es');end
    %OVERSAMPLING FOR ZOOMED ACQUISITIONS
    if Param.IsParameter('UGN1_GEO_sense_p_os_factor');Par.Mine.OverP=Param.GetValue('UGN1_GEO_sense_p_os_factor');end
    %SEGE FLAG TO ASSUME SECOND ECHO HAS THE SAME STRUCTURE AS FIRST ECHO
    if Param.IsParameter('UGN1_ACQ_epi_sege');Par.Mine.SAFE=Param.GetValue('UGN1_ACQ_epi_sege');else Par.Mine.SAFE='no';end
    %RF PHASE INCREMENT
    if Param.IsParameter('UGN1_ACQ_var_phaseIncr');Par.Mine.PhaseIncrement=Param.GetValue('UGN1_ACQ_var_phaseIncr');else Par.Mine.PhaseIncrement=0;end
        
    useGoalC=0;
    if useGoalC    
        %DUMPING THE GOAL-C PARAMETERS
        ParStruct=struct(ParStruct.Par40);
        lev1Fields=fieldnames(ParStruct);
        for n=1:length(lev1Fields)
            curFieldLev1=ParStruct.(lev1Fields{n});
            if isobject(curFieldLev1)
                L=length(curFieldLev1);%From here on we allocate arrays of structures. This was not done when creating the basic RF structure before, so we may encounter problems for arrays                        
                %Par.GoalC.(lev1Fields{n})=struct;%This throwed an error in the assignment below. Not sure how to preallocate memory
                for l=1:L
                   curFieldLev1scalar=curFieldLev1(l);                
                   if isobject(curFieldLev1scalar);curFieldLev1scalar=struct(curFieldLev1scalar);end
                   %Here another attempt to preallocate memory that has not provide much benefit
                   %if l==1;Par.GoalC.(lev1Fields{n})(1:L)=curFieldLev1scalar;else Par.GoalC.(lev1Fields{n})(l)=curFieldLev1scalar;end                
                   Par.GoalC.(lev1Fields{n})(l)=curFieldLev1scalar;
                end         
    %             %This was an alternative but has only stored the first element of the array of structs, not sure why
    %             lev2Fields=fieldnames(curFieldLev1);
    %             for s=1:length(lev2Fields)
    %                 Par.GoalC.(lev1Fields{n})(:).(lev2Fields{s})=curFieldLev1(:).(lev2Fields{s});                
    %             end            
            else
                Par.GoalC.(lev1Fields{n})=ParStruct.(lev1Fields{n});
            end
        end
    end
end
warning('on','MATLAB:structOnObject');



%This is an example:
%                     Rawfile: 'ba_29082017_1634145_2_2_pud2refneoheadV4.raw'
%               Parameter: [1×5491 GoalcParameter]
%                  Groups: [1×439 GoalcGroup]
%                 Objects: []
%                   Index: [1×22 struct]
%    ParAlphabeticalIndex: [1×15 struct]
%             NrParameter: 5491
%               NrObjects: 1092
%                NrGroups: 439
