function rec=InlineReconPipeline(twix_like)
% This is a pipeline modified for inline integration, it will be called by handle_connection_background_sense_acs

%RECONPIPELINE   Runs a reconstruction
%   RECONPIPELINE(CASEIN,BODYIN,REFEIN,FILEIN)
%   * CASEIN is the path with the data to reconstruct
%   * BODYIN is the file with the body coil data
%   * SURFIN is the file with the surface coil data
%   * DATAIN is the file with the actual data
%   ** REC is a recon structure
%

%ASSIGN NAMES
rec.Nam=struct('caseIn',twix_like.hdr.save_path,'bodyIn','','surfIn',twix_like.hdr.fileName,'dataIn',twix_like.hdr.fileName,...
'bodyInNoExt','','surfInNoExt',twix_like.hdr.fileName,'dataInNoExt',twix_like.hdr.fileName);

%ASSIGN PARAMETERS
rec=reconAlgorithm(rec);
typV=['B' 'S' 'x'];

    n0=2; % target scan with ACS, process ACS first as the reference, then the imaging data with SENSE recon

for n=n0:3
    %READ
    tsta=tic;  
    [rec,cont]=reconRead_inline(rec,n,twix_like);      
    tend=toc(tsta);
    fprintf('Time reading %s: %.2f\n',typV(n),tend);
    if cont;continue;end

    %INVERT
    tsta=tic;
    rec=reconInvert_inline(rec,n);
    fprintf('Time inverting %s: %.2f\n',typV(n),toc(tsta));
    
    %RECONSTRUCT
    tsta=tic;
    rec=reconReconstruct(rec,n);
    fprintf('Time reconstructing %s: %.2f\n',typV(n),toc(tsta));

    %WRITE
    tsta=tic;
    rec=reconWrite_inline(rec,n,'final');
    fprintf('Time writing %s: %.2f\n',typV(n),toc(tsta));
end



return
% % %%%%%HEREHEREHERE---CHECK HOW TO BUILD THE REFERENCE DATA
% % 
% % %BUILD THE REFERENCE DATA
% % if estSens
% %     for f=1:length(refIn{p})
% %         recS=[];
% %         fileName=strcat(pathIn{p},filesep,refIn{p}{f});
% %         TWS=mapVBVD(fileName);
% % 
% % 
% %         %%%%HEREHEREHERE---TO BE REMOVED, TO TEST GIULIO'S
% %         %rec.y=dynInd(TWS{2}.image,{':'},[16]);
% %         %x=dynInd(rec.y,2:2:size(rec.y,1),1);
% %         %x=permute(x,[1 3 4 2]);            
% %         %%x=gpuArray(x);
% %         %z=squeeze(multDimSum(abs(x).^2,[1 4]));
% %         %%z=gather(z);
% %         %zz=zeros(1,length(TWS{2}.image.Par));
% %         %for n=1:length(zz)
% %         %    zz(n)=z(TWS{2}.image.Lin(n),TWS{2}.image.Par(n));
% %         %end
% %         %size(zz)
% %         %figure
% %         %plot(log(zz))
% %         %%rec.Assign.z{2}=TWS{2}.image.Par;
% %         %%rec.Assign.z{3}=TWS{2}.image.Lin;
% % 
% % 
% %         recS=invert7T(TWS,supportReadS);
% %         recS=solveSensit7T(recS);
% %         save(strcat(fileName,'recS.mat'),'recS','-v7.3');
% %     end
% % end  
