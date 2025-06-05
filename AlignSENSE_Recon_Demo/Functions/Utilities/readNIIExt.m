function [x,MS,MT]=readNIIExt(dat,suff,gpu, underScore, applyRange)

%READNIIEXT   Is the readNII.m function to read nii files with extensions to allow
%suffixes that don't contain the underscore as well as using the NIFTI intensity rescaling option.
%   [X,MS,MT]=READNIIEXT(DAT,SUFF,{GPU},{UNDERSCORE},{APPLYRANGE})
%   * DAT indicates the path and name of the file
%   * SUFF is a list of suffixes of the files to be read
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu (0) arrays (empty, default depending on machine)
%   * {UNDERSCORE} is a flag to include an underscaore to the suffixes (Defaults to 0).
%   * {APPLYRANGE} is a flag to use the NIFTI intensity rescaling option (Defaults to 1).
%   ** X is the returned data
%   ** MS is the returned spacing
%   ** MT is the returned orientation
%
%   Yannick Brackenier 2022-01-30

if nargin<2 || isempty(suff);suff={''};end
if nargin<3 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end
if nargin<4 || isempty(underScore);underScore=0;end
if nargin<5 || isempty(applyRange);applyRange=1;end

N=length(suff);
x=cell(1,N);MS=cell(1,N);MT=cell(1,N);
for n=1:N
    if underScore%YB modification
        fileAbs{1}=sprintf('%s_%s.nii',dat,suff{n});
        fileAbs{2}=sprintf('%s_%s.nii.gz',dat,suff{n});
    else
        fileAbs{1}=sprintf('%s%s.nii',dat,suff{n});
        fileAbs{2}=sprintf('%s%s.nii.gz',dat,suff{n});
    end
    
    for l=1:length(fileAbs)
        if exist(fileAbs{l},'file');nii=load_untouch_nii(fileAbs{l});break;end
        if l==length(fileAbs);error('File %s does not exist',fileAbs{1});end
    end
    if l==1;filePhase=sprintf('%s%sPh.nii',dat,suff{n});else filePhase=sprintf('%s%sPh.nii.gz',dat,suff{n});end
                
    x{n}=nii.img;
    
    if applyRange %YB modification
        if isfield(nii.hdr.dime,'scl_slope') && ~isempty(nii.hdr.dime.scl_slope) && nii.hdr.dime.scl_slope>0; x{n} = x{n} * nii.hdr.dime.scl_slope;end
        if isfield(nii.hdr.dime,'scl_inter') && ~isempty(nii.hdr.dime.scl_inter) && nii.hdr.dime.scl_inter~=0; x{n} = x{n} + nii.hdr.dime.scl_inter;end
    end
    if gpu;x{n}=gpuArray(x{n});end
    if exist(filePhase,'file')
        nii=load_untouch_nii(filePhase);
        p=nii.img;
        if gpu;p=gpuArray(p);end
        x{n}=x{n}.*exp(1i*p);
    end
    MS{n}=nii.hdr.dime.pixdim(2:4);
    MT{n}=eye(4);MT{n}(1,:)=nii.hdr.hist.srow_x;MT{n}(2,:)=nii.hdr.hist.srow_y;MT{n}(3,:)=nii.hdr.hist.srow_z;
    %YB: check if empty s_form
    %if all(MT{n}(:)==0); MT{n} = convertNIIGeom(MS{n}, [nii.hdr.hist.quatern_b,nii.hdr.hist.quatern_c,nii.hdr.hist.quatern_d ], 'qForm', 'sForm');end

end
