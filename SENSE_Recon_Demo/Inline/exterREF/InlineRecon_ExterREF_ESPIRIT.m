function rec=InlineRecon_ExterREF_ESPIRIT(twix_like)
% This is a pipeline modified for inline integration, it will be called by handle_connection_background_sense_acs

%RECONPIPELINE   Runs a reconstruction
%   RECONPIPELINE(CASEIN,BODYIN,REFEIN,FILEIN)
%   * CASEIN is the path with the data to reconstruct
%   * BODYIN is the file with the body coil data
%   * SURFIN is the file with the surface coil data
%   * DATAIN is the file with the actual data
%   ** REC is a recon structure
%

% %ASSIGN NAMES
rec.Nam=struct('caseIn',twix_like.hdr.save_path,'bodyIn','','surfIn',twix_like.hdr.fileName,'dataIn',twix_like.hdr.fileName,...
'bodyInNoExt','','surfInNoExt',twix_like.hdr.fileName,'dataInNoExt',twix_like.hdr.fileName);

%ASSIGN PARAMETERS
rec=reconAlgorithm(rec);
typV=['B' 'S' 'x'];

% if isempty(rec.Nam.bodyIn)
n=2; % target scan with ACS, process ACS first as the reference, then the imaging data with SENSE recon
%     fprintf('No body coil\n');
% end

%READ
rec.Alg.useBuiltInCalibration = -2; % calculate the csm by the external REF
tsta=tic;  
[rec,cont]=reconRead_inline(rec,n,twix_like);     
tend=toc(tsta);
fprintf('Time reading %s: %.2f\n',typV(n),tend);


return
