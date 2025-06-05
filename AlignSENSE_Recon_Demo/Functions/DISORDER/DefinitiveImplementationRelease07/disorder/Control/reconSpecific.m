function rec=reconSpecific(rec)

%RECONSPECIFIC   Overwrites Plan and Alg parameters for specific purposes
%   REC=RECONSPECIFIC(REC)
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (.Names), the status of the reconstruction (.Fail) 
%   and the .lab information (.Par)
%   ** REC is a reconstruction structure with filled plan information 
%   (.Plan for fixed plan information and .Dyn for dynamic plan 
%   information)
%

if rec.Fail;return;end

if isfield(rec.Names,'Specific') && ~isempty(rec.Names.Specific)
    if strcmp(rec.Names.Specific,'Detuning')
        rec.Alg.DetectAnomalies=7;rec.Dyn.Debug=0;
    end
    if strcmp(rec.Names.Specific,'Spikes')
        rec.Alg.DetectAnomalies=4;rec.Dyn.Debug=0;
    end
    if strcmp(rec.Names.Specific,'WithinShot')
        rec.Plan.SuffOu='FullWithin';
        rec.Alg.AlignedRec=4;rec.Alg.parXT.NWend=16;rec.Alg.parXT.accel=[0 0];
        rec.Alg.parXT.saveFinal=1;
    end
    if strcmp(rec.Names.Specific,'WriteDISORDER')
        rec.Alg.parXT.saveFinal=1;
    end
end

%This was before in ReconPlanner, but it is better run after potential
%modifications
if rec.Dyn.Debug>=2
    fprintf('Expected size of reconstructed grid:%s\n',sprintf(' %d',rec.Dyn.ExpSiz));
    fprintf('Use of GPU acceleration: %d\n',rec.Dyn.GPU);
end        

