function reference = referenceFromReconData_bucket(recon_data, connection)
    
    hdr = recon_data.bits.buffer.bucket_headers; % ZN: small headers of every kspace line
    samplingInfo = recon_data.bits.buffer.sampling_description;
    
    % ZN: to extract the right header (of the proper lines) to generate the
    % reference -> header of the reconstructed image
    % this function works if your readin is bucket data
    % refer to referenceFromReconData function if your readin is buffer
    %
    % Here, you need to find the first sampled data in gridded kspace, and
    % take it small header to generate reference
    timeStamp = ( squeeze(hdr.acquisition_time_stamp) );
    idxSamples = find(timeStamp >0 );
    kspace_1_extracted = hdr.kspace_encode_step_1(idxSamples);
    kspace_2_extracted = hdr.kspace_encode_step_2(idxSamples);
    min_kspace_1 = min(kspace_1_extracted);
    idxKsapce_1_extract = find(kspace_1_extracted == min_kspace_1);
    idxKspace_1 = idxSamples(idxKsapce_1_extract);
    kspace_2_tmp = kspace_2_extracted(idxKsapce_1_extract);
    min_kspace_2 = min(kspace_2_extracted);
    idxFirstSample = idxKspace_1(kspace_2_tmp == min_kspace_2);
    if ~isempty(idxFirstSample)
        reference = structfun(@(field) field(:, idxFirstSample(1)), hdr, 'UniformOutput', false);
    else
        fprintf('Cannot find the idx of the first sampled data in gridded kspace somehow \n');
        fprintf('Take the first sampled data in the timestamp for reference generation \n');
        fprintf('!!!! Might be error to retrieve reconstructed images to the scanner \n');
        reference = structfun(@(field) field(:, idxSamples(1)), hdr, 'UniformOutput', false);
    end
    
   
    %%% ad fov since not in hdr
    rec_FOVy = single(samplingInfo.recon_fov(2));
    rec_FOVz = single(samplingInfo.recon_fov(3));
    rec_FOVx = single(samplingInfo.recon_fov(1));
    reference.field_of_view = [rec_FOVx rec_FOVy rec_FOVz];  %FOVRec
        
    fprintf('--- referenceFromReconData:: need to convert orientation properly.')
end