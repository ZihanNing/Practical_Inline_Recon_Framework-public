function x_gibbsringredu = gibbs_ringing_reduction_mrtrix3(outDir, filename)
    %%%%%%%%%%%%%%%%%%%%%%
    % THIS FUNCTION WILL CALL MRTRIX3 FOR GIBBS RINGING REDUCTION
    % PLEASE MAKE SURE THAT THE MRTRIX3 HAS BEEN INSTALLED PROPERLLY
    % OUTDIR: DIR OF THE TO-BE-PROCESSED NIFTI FILE
    % FILENAME: FILENAME OF THE TO-BE-PROCESSED IMAGE FILE (NIFTI)
    % X_GIBBSRINGREDU: THE IMAGE AFTER GIBBS RINGING REDUCTION
    % THE PROCESSED NIFTI FILE WILL BE SAVED IN THE SAVE FOLDER (_gibbsringredu)
    %
    % BY ZIHAN N
    % 11-SEP-2024
    %%%%%%%%%%%%%%%%%%%%%%

    % Define the bash script
    bashScript = 'Bash_script/run_mrdegibbs.sh';
    
    % Build the command to call the bash script with the directory and filename
    command = sprintf('./%s %s %s', bashScript, outDir, filename);
    
    % Run the bash script and capture the status
    [status, cmdOut] = system(command);
    
    % Check the returned status from the bash script
    if status == 0
        % Success: Load the generated NIfTI file into MATLAB
        generatedFile = fullfile(outDir, [filename '_gibbsringredu.nii']);
        disp('File processed successfully. Loading NIfTI file...');
        x_gibbsringredu = niftiread(generatedFile);
        
        % Display or further process the NIfTI data if needed
        disp('NIfTI file loaded into MATLAB:');
    else
        % Failure: Print the error message
        disp('Bash script failed to process the image:');
        disp(cmdOut);
    end
end
