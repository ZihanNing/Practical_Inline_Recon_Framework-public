function hdr_image = referenceFromReconData_general(twix)
%%% This function generate the header structure for reconstructed image by
%%% the custom recon algorithm by the twix-like raw header (convereted from
%%% ISMRMRD raw)
%%%
%%% INPUT:
%%% twix: the twix-like structure (converted from the ISMRMRD raw), you can
%%% remove the dataset in this strcture (twix.data,twix.reference,twix.noise, etc)
%%% But make sure that your twix-like structure has the following header
%%% info: twix.hdr,twix.sampling_description,twix.bucket_header
%%% 
%%% OUTPUT:
%%% hdr_image: the header of the reconstructed image
%%%
%%% by Zihan @ king's
%%% March-2025
    
    hdr = twix.bucket_header; % ZN: small headers of every kspace line
    samplingInfo = twix.sampling_description;
    
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
        hdr_image = structfun(@(field) field(:, idxFirstSample(1)), hdr, 'UniformOutput', false);
    else
        fprintf('Cannot find the idx of the first sampled data in gridded kspace somehow \n');
        fprintf('Take the first sampled data in the timestamp for reference generation \n');
        fprintf('!!!! Might be error to retrieve reconstructed images to the scanner \n');
        hdr_image = structfun(@(field) field(:, idxSamples(1)), hdr, 'UniformOutput', false);
    end
    
   
    %%% ad fov since not in hdr
    rec_FOVy = single(samplingInfo.recon_fov(2));
    rec_FOVz = single(samplingInfo.recon_fov(3));
    rec_FOVx = single(samplingInfo.recon_fov(1));
    hdr_image.field_of_view = [rec_FOVx rec_FOVy rec_FOVz];  %FOVRec
end