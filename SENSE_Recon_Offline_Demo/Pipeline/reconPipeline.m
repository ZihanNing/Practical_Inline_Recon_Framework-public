function rec=reconPipeline(caseIn,bodyIn,surfIn,dataIn)

%RECONPIPELINE   Runs a reconstruction
%   RECONPIPELINE(CASEIN,BODYIN,REFEIN,FILEIN)
%   * CASEIN is the path with the data to reconstruct
%   * BODYIN is the file with the body coil data
%   * SURFIN is the file with the surface coil data
%   * DATAIN is the file with the actual data
%   ** REC is a recon structure
%

%ASSIGN NAMES
rec.Nam=struct('caseIn',caseIn,'bodyIn',bodyIn,'surfIn',surfIn,'dataIn',dataIn,...
'bodyInNoExt',removeExtension(bodyIn,'.dat'),'surfInNoExt',removeExtension(surfIn,'.dat'),'dataInNoExt',removeExtension(dataIn,'.dat'));

%ASSIGN PARAMETERS
rec=reconAlgorithm(rec);
typV=['B' 'S' 'x'];

if isempty(rec.Nam.bodyIn)
    n0=2;
    fprintf('No body coil\n');
end

for n=n0:3
    if n==1 && exist(fullfile(caseIn,sprintf('%s_B.nii.gz',rec.Nam.bodyInNoExt)),'file') && rec.Pip.readFromFile(1)
        [rec.B,rec.Geom.B.MS,rec.Geom.B.MT]=readnii(fullfile(caseIn,sprintf('%s_B',rec.Nam.bodyInNoExt)));
        fprintf('Reading %s from file\n',typV(n));
        continue;
    end
    if n==2 && exist(fullfile(caseIn,sprintf('%s_S.nii.gz',rec.Nam.surfInNoExt)),'file') && rec.Pip.readFromFile(2)
        [rec.S,rec.Geom.S.MS,rec.Geom.S.MT]=readnii(fullfile(caseIn,sprintf('%s_S',rec.Nam.surfInNoExt)));       
        rec.M=readnii(fullfile(caseIn,sprintf('%s_M',rec.Nam.surfInNoExt)));
        rec.B=readnii(fullfile(caseIn,sprintf('%s_B',rec.Nam.surfInNoExt)));
        fprintf('Reading %s from file\n',typV(n));
        continue;
    end
    if n==3 && exist(fullfile(caseIn,sprintf('%s_S.nii.gz',rec.Nam.dataInNoExt)),'file') && rec.Pip.readFromFile(3) && rec.Alg.useBuiltInCalibration
        if rec.Alg.useBuiltInCalibration==2
            [rec.S,rec.Geom.S.MS,rec.Geom.S.MT]=readnii(fullfile(caseIn,sprintf('%s_S',rec.Nam.dataInNoExt)));
            rec.M=readnii(fullfile(caseIn,sprintf('%s_M',rec.Nam.dataInNoExt)));
            rec.B=readnii(fullfile(caseIn,sprintf('%s_B',rec.Nam.dataInNoExt)));
            rec.xS=readnii(fullfile(caseIn,sprintf('%s_xS',rec.Nam.dataInNoExt)));
            NC=size(rec.B,4);NS=size(rec.S);
            rec.S=reshape(rec.S,[NS(1:3) NS(4)/NC NC]);
        end
        if exist(fullfile(caseIn,sprintf('%s_PS.nii.gz',rec.Nam.dataInNoExt)),'file');[rec.PS,rec.Geom.PS.MS,rec.Geom.PS.MT]=readnii(fullfile(caseIn,sprintf('%s_PS',rec.Nam.dataInNoExt)));end
        fprintf('Reading %s from file\n',typV(n));
    end
    fprintf('Processing %s\n',typV(n));

    %READ
    tsta=tic;  
    [rec,cont]=reconRead(rec,n);    
    tend=toc(tsta);
    fprintf('Time reading %s: %.2f\n',typV(n),tend);
    if cont;continue;end

    %INVERT
    tsta=tic;
    rec=reconInvert(rec,n);
    fprintf('Time inverting %s: %.2f\n',typV(n),toc(tsta));
    
    %RECONSTRUCT
    tsta=tic;
    rec=reconReconstruct(rec,n);
    fprintf('Time reconstructing %s: %.2f\n',typV(n),toc(tsta));

    %WRITE
    tsta=tic;
    rec=reconWrite(rec,n);
    fprintf('Time writing %s: %.2f\n',typV(n),toc(tsta));
end

if nargout<1;rec=[];end

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