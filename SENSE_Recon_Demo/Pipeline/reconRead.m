function [rec,cont]=reconRead(rec,typ)

%RECCONREAD   Reads data
%   REC=RECONREAD(REC,TYP)
%   * REC is a recon structure
%   * TYP is the type of data to read, one of the following, 1->body, 
%   2->surface, 3->data
%   ** REC is a recon structure
%

if typ==1;fileName=rec.Nam.bodyIn;%Body
elseif typ==2;fileName=rec.Nam.surfIn;%Surface
elseif typ==3;fileName=rec.Nam.dataIn;%Data
end

ND=16;
typV=['B' 'S' 'x'];
t=typV(typ);
fileName=fullfile(rec.Nam.caseIn,filesep,fileName);

%We use previously read data if surface coil and actual data are the same
%if typ==3 && strcmp(fileName,fullfile(rec.Nam.caseIn,filesep,rec.Nam.surfIn))
%    rec.TW.(t)=rec.TW.S;
%    rec.TWD.(t)=rec.TWD.S;
%else    
    rec.TW.(t)=mapVBVD(fileName);
    if iscell(rec.TW.(t));
        fprintf('Cell information, taking last element\n');
        for m=1:length(rec.TW.(t));fprintf('Number of profiles element %d: %d\n',m,rec.TW.(t){m}.image.NAcq);end
        rec.TW.(t)=rec.TW.(t){length(rec.TW.(t))};        
    end%We take the last element if it is a cell
    rec.TW.(t).image.flagIgnoreSeg=1;%Ignores segment mdh index
    %return
        
    if isfield(rec.TW.(t),'refscan') && ~isfield(rec,'PS') && rec.Alg.useBuiltInCalibration>0 && ~strcmp(t,'x')
        %RENAME FIELD AND READ DATA
        if rec.Alg.useBuiltInCalibration==2
            typ=2;
            rec.TW=renameStructField(rec.TW,'x','S');t=typV(typ);
            rec=rmfield(rec,'B');
        end
        if rec.Alg.removeOversamplingReadout;rec.TW.(t).refscan.flagRemoveOS=true;end
        rec.TWD.(t)=dynInd(rec.TW.(t).refscan,{':'},ND);      
        %size(rec.TWD.(t))
    
        fprintf('----------------------\n');
        fprintf('Reconstructing refscan\n');tsta=tic;

        %INVERT
        tsta=tic;
        rec=reconInvert(rec,typ,'refscan');
        fprintf('Time inverting %s: %.2f\n',t,toc(tsta));
    
        %RECONSTRUCT
        tsta=tic;
        rec=reconReconstruct(rec,typ,'refscan');
        fprintf('Time reconstructing %s: %.2f\n',t,toc(tsta));


        fprintf('----------------------\n');        

        if rec.Alg.useBuiltInCalibration==2
            typ=3;
            rec.TW=renameStructField(rec.TW,'S','x');t=typV(typ);
        end
        cont=1;
    else
        cont=0;
    end

    %READ DATA
    %rec.TWD.(t)=dynInd(rec.TW.(t).image,{':'},ND);
    %rec.TWD.(t)=dynInd(rec.TW.(t).image,{':'},7);
    if rec.Alg.removeOversamplingReadout;rec.TW.(t).image.flagRemoveOS = true;end
    %rec.TWD.(t)=rec.TW.(t).image(:,:,:,:,:,:,:,1);
    rec.TWD.(t)=dynInd(rec.TW.(t).image,{':'},ND);
%end
if typ>1 && isfield(rec,'TWD') && isfield(rec.TWD,typV(typ-1));rec.TWD=rmfield(rec.TWD,typV(typ-1));end%We release data
