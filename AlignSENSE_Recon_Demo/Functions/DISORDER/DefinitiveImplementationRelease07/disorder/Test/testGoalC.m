%load /home/lcg13/Data/tu_01062019_1134218_8_2_mv3ssfptr6fa46rf292senseV4.mat
%Par=Par.GoalC.Parameter;
%load /home/lcg13/Data/tu_01062019_1137215_9_2_mv3ssfptr6fa12rf95senseV4.mat
%Par=Par.GoalC.Parameter;
load /home/lcg13/Data/tu_01062019_1147058_12_2_mv3spgrtr8p2fa12rf50senseV4.mat
Par=Par.GoalC.Parameter;


for n=1:length(Par)
    if ~iscell(Par(n).values(1)) && Par(n).values(1)==50% && (strcmp(Par(n).name,'UGN1_ACQ_var_phaseIncr') || strcmp(Par(n).name,'MP_FFE_ssfp_phase_step') || strcmp(Par(n).name,'EX_ACQ_var_phaseIncr'))
        Par(n).name
        Par(n).values
    end
end