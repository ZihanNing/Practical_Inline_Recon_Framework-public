function rec=reconWrite_inline(rec,typ,status)

%RECCONWRITE   Reads data
%   RECONWRITE(REC,TYP,{FIELD})
%   * REC is a recon structure
%   * TYP is the type of data to read, one of the following, 1->body, 
%   2->surface, 3->data
%   * {FIELD} is the field with header info, it defaults to 'image', other
%   options are 'refscan'
%
%   * {Status} the current status of calculation
%       - reconInvert: called after reconInvert


if typ==1
    fileName=rec.Nam.bodyInNoExt;%Body
else % ZN: refscan or high-res
    fileName=rec.Nam.dataInNoExt;%Data
end

typV=['B' 'S' 'x'];
t=typV(typ);
fileName=fullfile(rec.Nam.caseIn,fileName);

if typ == 2 % ZN: refscan
    MSS=rec.Geom.S.MS;
    MTT=rec.Geom.S.MT;
else % ZN: high-res
    MT=rec.Geom.(t).MT;
    MS=rec.Geom.(t).MS;
end

if typ==2 % ZN: refscan, to be modified
    switch status
        case 'reconInvert'
            dims = ndims(rec.(t)); 
            if dims>4
                idx = repmat({':'},1,dims);
                for k = 5:dims; idx{k} = 1;end % ZN: took first echo/inv
                save_S = RSOS(rec.(t)(idx{:}));
            else
                save_S = RSOS(rec.(t));
            end
            writenii(fileName,save_S,'S_reconInvert',MSS,MTT); % ZN: write the refscan image after reconInvert 
            fprintf('nii for high-res dataset after reconInvert (aliased) saved \n');
        case 'refscan'
            writenii(fileName,rec.M,'M',MSS,MTT);
            writenii(fileName,rec.B,'B',MSS,MTT);
            writenii(fileName,rec.S,'S',MSS,MTT);
    end
else % ZN: high-res, to be modified
    switch status
        case 'reconInvert'
            save_x = ones([size(rec.x,1:3),size(rec.x,8)]); % save all echoes
            for i = 1:size(rec.x,8)
                save_x(:,:,:,i) = RSOS(squeeze(rec.x(:,:,:,:,1,1,1,i)));
            end
            writenii(fileName,save_x,'x_reconInvert',MS,MT); % ZN: write the high-res image before SENSE recon (aliased) 
            fprintf('nii for high-res dataset after reconInvert (aliased) saved \n');
        case 'final'
            writenii(fileName,squeeze(rec.x),'x',MS,MT);
    end
end

