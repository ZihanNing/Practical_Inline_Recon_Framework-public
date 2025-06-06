
function TW = dat2TWIX(fileName)
     
%   Yannick Brackenier

%%% NAME HANDLING
[path,file,~] = fileparts(fileName);
fileName=fullfile(path,file);%Remove the suffix in case it does not have the .dat extension

%%% CONVERT
if existsFileVar(strcat(fileName,'.dat'))~=1
    TW=[];
    warning('dat2TWIX: %s does not exist. Empty twix object returned.', strcat(fileName,'.dat'));
else
    evalc('TW = mapVBVD(strcat(fileName,''.dat''))');%Suppress command line output
end

