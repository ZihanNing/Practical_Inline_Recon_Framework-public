function lock_check()

	lock='/tmp/pause_recon';
	acknowledge='/tmp/recon_paused';        

    c=0;
	while exist(lock, 'file') == 2
		if exist(acknowledge, 'file') ~= 2
            if gpuDeviceCount>0;dev=gpuDevice;else dev=[];end
            reset(dev);
			disp(['creating ' acknowledge])
			fclose(fopen(acknowledge, 'w'));
			fileattrib(acknowledge,'+w','a')
        end
        if mod(c,360)==0;disp(['waiting for lock file to be deleted: ' lock]);end%This is written roughly once per hour
        c=c+1;
		pause(10)
	end

	if exist(acknowledge, 'file') == 2
		disp(['deleting ' acknowledge])
	  	delete(acknowledge);
	end

end