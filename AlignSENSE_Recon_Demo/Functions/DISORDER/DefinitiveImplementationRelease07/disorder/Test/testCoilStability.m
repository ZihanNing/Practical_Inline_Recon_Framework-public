inpIm{1}='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Coils/Adult/2019_01_28/Ad_76430/Re-Se/ad_28012019_1136352_3_2_dhcp8refneoheadV4';
inpIm{2}='/home/lcg13/Data/pnrawDe/ReconstructionsDebug06/Coils/Adult/2019_01_28/Ad_76430/Re-Se/ad_28012019_1250407_21_2_dhcp8refhead32V4';
suff{1}='Se';
gpu=1;
for n=1:2
    [y,MS,MT]=readNII(inpIm{n},suff,gpu);
    x{n}=y{1};
end

%FAILS BECAUSE NOT SAME RESOLUTION
for n=1:2
    z{n}=normm(x{2*n},x{2*n-1},4)./normm(x{2*n},[],4);
    %writeNII(inpIm{2*n},{'CoilDifference'},{z{n}});
end

visReconstruction(z{1})