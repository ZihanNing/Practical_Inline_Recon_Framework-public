function rec=reconWrite(rec,typ,field)

%RECCONWRITE   Reads data
%   RECONWRITE(REC,TYP,{FIELD})
%   * REC is a recon structure
%   * TYP is the type of data to read, one of the following, 1->body, 
%   2->surface, 3->data
%   * {FIELD} is the field with header info, it defaults to 'image', other
%   options are 'refscan'
%

if nargin<3 || isempty(field);field='image';end

if typ==1
    fileName=rec.Nam.bodyInNoExt;%Body
elseif typ==2
    if strcmp(field,'image');fileName=rec.Nam.surfInNoExt;%Surface
    elseif strcmp(field,'refscan');fileName=rec.Nam.dataInNoExt;
    else error('Unrecognized field: %s',field);
    end
elseif typ==3;fileName=rec.Nam.dataInNoExt;%Data
end

typV=['B' 'S' 'x'];
t=typV(typ);
fileName=fullfile(rec.Nam.caseIn,fileName);

MT=rec.Geom.(t).MT;
MS=rec.Geom.(t).MS;
if strcmp(field,'refscan')
    MSS=rec.Geom.S.MS;
    MTT=rec.Geom.S.MT;
end

if typ<3 || strcmp(field,'refscan')
    if ~strcmp(field,'refscan');writenii(fileName,rec.(t),t,MS,MT);
    else writenii(fileName,rec.S,'S',MSS,MTT);
    end
    if typ==2 || strcmp(field,'refscan')
        if strcmp(field,'refscan');writenii(fileName,rec.M,'M',MSS,MTT);else writenii(fileName,rec.M,'M',MS,MT);end
        if strcmp(field,'refscan');writenii(fileName,rec.B,'B',MSS,MTT);else writenii(fileName,rec.B,'B',MS,MT);end
        %return
        if isfield(rec,'W') && ~isempty(rec.W)
            if strcmp(field,'refscan');writenii(fileName,rec.W,'W',MSS,MTT);else writenii(fileName,rec.W,'W',MS,MT);end
        end
        writenii(fileName,rec.xS,'xS',MS,MT);
        if isfield(rec,'PS');writenii(fileName,rec.PS,'PS',MS,MT);end
        %rec=rmfield(rec,'xS');
    end
else
    if rec.Alg.writeRaw;writenii(fileName,x,'x',MS,MT,gatherStruct(rec));
    else 
        NX=size(rec.x);
        rec.x=resampling(rec.x,[NX(1) rec.Enc.x.OutN(2:3)],3);
        tI=MT*[0;-0.5*(rec.Enc.x.OutN(2:3)-rec.Enc.x.RecN(2:3))';1];
        MT(1:3,4)=tI(1:3);
        writenii(fileName,rec.x,'x',MS,MT);
    end
end
