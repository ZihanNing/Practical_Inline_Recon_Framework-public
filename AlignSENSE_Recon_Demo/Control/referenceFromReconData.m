
function reference = referenceFromReconData(recon_data, connection)
    
    hdr = recon_data.bits.buffer.headers;
    samplingInfo = recon_data.bits.buffer.sampling_description;
    
    timeStamp = ( squeeze(hdr.acquisition_time_stamp) );
    idxSamples = find(timeStamp >0 );
    idxValidSample = idxSamples(1);
    
    reference = structfun(@(field) field(:, idxValidSample), hdr, 'UniformOutput', false);
    
   
    %%% ad fov since not in hdr
    %FOVRec
    rec_FOVx = single(samplingInfo.recon_fov(1));
    rec_FOVy = single(samplingInfo.recon_fov(2));
    rec_FOVz = single(samplingInfo.recon_fov(3));
    reference.field_of_view = [rec_FOVx rec_FOVy rec_FOVz];
    
    fprintf('--- referenceFromReconData:: need to convert orientation properly.')
end



