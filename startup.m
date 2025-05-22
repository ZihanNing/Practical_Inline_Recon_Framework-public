cd('/home/zn23/gadgetron_test')
addpath(genpath('/home/zn23/MATLAB Add-Ons/Toolboxes'))

warning off parallel:gpu:device:DeviceLibsNeedsRecompiling
try
    gpuArray.eye(2)^2;
catch ME
end
% try
%     nnet.internal.cnngpu.reluForward(1);
% catch ME
% end
