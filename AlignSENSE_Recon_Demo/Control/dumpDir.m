

function [dirName] = dumpDir()

% dirName = fullfile(cd,'DATA_DUMP'); % ZN: creat DATA_DUMP in the DISORDER_gadgetron folder
dirName = '/home/gadgetron/data'; % ZN: relocate the DATA_DUMP to another disk
if ~isfolder(dirName); mkdir(dirName);end