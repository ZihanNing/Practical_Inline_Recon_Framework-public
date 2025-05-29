addpath(genpath('.'))

if useGPU;gpuDevice(2);end

idData=1;%[2 5];%[1 3 4];%[];%1:3;%[];%Reconstruction to run, if empty all

%FOLDER AND FILES TO RUN
siemensStudies;
if ~isempty(idData)
    dataIn=dataIn(idData);
    surfIn=surfIn(idData);
end

%dataIn=[];

%RECONSTRUCTION PIPELINE
for n=1:length(dataIn)
    rec=[];
    rec=reconPipeline(caseIn,bodyIn,surfIn{n},dataIn{n});
%     %continue
%     if strcmp(dataIn{n},'meas_MID00751_FID57555_SWI_LinearCarte_separateACS.dat') || ... %Adjustment of FOV
%        strcmp(dataIn{n},'meas_MID00750_FID57554_SWI_DISORDER_still.dat') || ...
%        strcmp(dataIn{n},'meas_MID00743_FID57547_SWI_DISORDER_move.dat') || ...
%        strcmp(dataIn{n},'meas_MID00107_FID16038_T1_MPRAGE_DISORDER_twinsACS.dat') 
%        %strcmp(dataIn{n},'meas_MID00742_FID57546_SWI_LinearCarte_integratedACS.dat') || ...
%         fileName=rec.Nam.dataInNoExt;%Data
%         fileName=fullfile(rec.Nam.caseIn,fileName);
%         fileNameAux=strcat(fileName,'_x.nii.gz');
%         if isfile(fileNameAux)
%             [x,MS,MT]=readnii(fileNameAux);
%             if strcmp(dataIn{n},'meas_MID00751_FID57555_SWI_LinearCarte_separateACS.dat')
%                 [xAux,MSAux,MTAux]=readnii(fullfile(rec.Nam.caseIn,'nifti_results','case1','SWI_Lin_still_offline_separateACS_results','meas_MID00751_FID57555_SWI_UKBB_separateREF_still_Aq_MotCorr_echo1'));
%             %elseif strcmp(dataIn{n},'meas_MID00742_FID57546_SWI_LinearCarte_integratedACS.dat')
%             %    [xAux,MSAux,MTAux]=readnii(fullfile(rec.Nam.caseIn,'nifti_results','case1','SWI_LinearCarte_integratedACS','matlab_SENSE','meas_MID00428_FID19659_SWI_UKBB_Aq_MotCorr_echo1.nii'));            
%             elseif strcmp(dataIn{n},'meas_MID00750_FID57554_SWI_DISORDER_still.dat')
%                 [xAux,MSAux,MTAux]=readnii(fullfile(rec.Nam.caseIn,'nifti_results','case1','SWI_DISORDER_MoCo_still_results','meas_MID00750_FID57554_SWI_DISORDER_R2_TS58_still_Aq_MotCorr_echo1'));
%             elseif strcmp(dataIn{n},'meas_MID00743_FID57547_SWI_DISORDER_move.dat')
%                 [xAux,MSAux,MTAux]=readnii(fullfile(rec.Nam.caseIn,'nifti_results','case1','SWI_DISORDER_MoCo_move_results','meas_MID00743_FID57547_SWI_DISORDER_R2_TS58_Gad_matched_Aq_MotCorr_echo1'));
%             elseif strcmp(dataIn{n},'meas_MID00107_FID16038_T1_MPRAGE_DISORDER_twinsACS.dat')
%                 [xAux,MSAux,MTAux]=readnii(fullfile(rec.Nam.caseIn,'nifti_results','MPRAGE','meas_MID00107_FID16038_T1_MPRAGE_DISORDER_twinsACS_Aq_MotCorr'));
%             end
%             MTNew=MT;
%             MTNew(1:3,4)=MTAux(1:3,4);
%             fileNameAux=strcat(fileName,'_xMapped.nii.gz');
%             writenii(fileNameAux,x,[],MS,MTNew);                
%         end
%     end
end

return

