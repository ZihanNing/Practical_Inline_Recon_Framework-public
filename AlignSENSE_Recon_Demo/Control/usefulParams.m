
%bufferData.is_flag_set(bufferData.ACQ_IS_NOISE_MEASUREMENT)

%noise_dwell_time = acquisition.header.sample_time_us;

%to increase series of images in .h5
%images(1).header.image_series_index=images(1).header.image_series_index+1;

%if acquisitions(counter).is_flag_set(acquisitions(1).ACQ_LAST_IN_MEASUREMENT)
%	break; 
%end