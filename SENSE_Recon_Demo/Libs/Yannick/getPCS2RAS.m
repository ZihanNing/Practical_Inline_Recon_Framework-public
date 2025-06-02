

function [PCS2RAS] = getPCS2RAS ( patientOrientation )

%TO COMPLETE: see Siemens IDEA manual

if strcmp(patientOrientation,'HFS')
    PCS2RAS =  diag([-1 -1 +1 1]);
elseif strcmp(patientOrientation,'FFS')
    PCS2RAS =  diag([+1 -1 -1 1]);
else
    error('getPCS2RAS:: patientOrientation not implemented yet. Time to look it up now in the manual!!')
end