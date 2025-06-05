function changePermissions(path,rootOu)

%CHANGEPERMISSIONS   Changes permissions of specific output folder
%   CHANGEPERMISSIONS(PATH,ROOTOU)
%   * PATH is the relative path to raw data
%   * {ROOTOU} is the root of the destination folder for the parsed data 
%   (and later for the reconstructions).
%

[user,versCode]=versRecCode;
if nargin<2 || isempty(rootOu);rootOu=fullfile(filesep,'home',user,'Data','rawDestin',sprintf('Reconstructions%s',versCode));end

day=strsplit(path,'/');day=day{1};day=fullfile(rootOu,day);
systCall=sprintf('chgrp -R pnraw_dhcp_recon_releases %s',day);
system(systCall);
systCall=sprintf('chmod -R 775 %s',day);
system(systCall);
systCall=sprintf('chmod -R -x+X %s',day);
system(systCall);
