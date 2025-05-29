function [x,MS,MT,datima,hdrU]=readnii(fil,suff,gpu,headonly)

%READNII   Reads a set of nii files and potentially associated matlab files
%   [X,MS,MT,DATIMA,HDRU,HDRT]=READNII(FIL,{SUFF},{GPU},{HEADONLY})
%   * FIL indicates the path and name of the file
%   * {SUFF} is a list of suffixes of the files to be read, if empty no
%   suffix is read and result is returned as an array instead of a cell
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu 
%   (0) arrays (empty, default depending on machine)
%   * {HEADONLY} allows to read header only, defaults to 0
%   ** X is the image data
%   ** MS is the spacing
%   ** MT is the orientation
%   ** DATIMA is the associated metadata
%   ** HDRU is the nifti header (untouched)
%

if nargin<2;suff=[];end
if nargin<3 || isempty(gpu);gpu=useGPU;end
if nargin<4 || isempty(headonly);headonly=0;end

%IF EXTENSION IS ENTERED, WE REMOVE IT
fil=removeExtension(fil,{'.nii','.gz'});

N=length(suff);
Nc=max(N,1);
x=cell(1,Nc);MS=cell(1,Nc);MT=cell(1,Nc);datima=cell(1,Nc);hdrU=cell(1,Nc);
for n=1:Nc
    if N>0;suff{n}=strcat('_',suff{n});else suff{n}='';end
    fileAbs{1}=strcat(fil,suff{n},'.nii');
    fileAbs{2}=strcat(fil,suff{n},'.nii.gz');
    for l=1:length(fileAbs)
        if exist(fileAbs{l},'file')
            if headonly;nii.hdr=load_untouch_header_only(fileAbs{l});nii.img=[];else nii=load_untouch_nii(fileAbs{l});end
            if l==1;outExt='.nii';else outExt='.nii.gz';end
            break;
        end
        if l==length(fileAbs);error('File %s does not exist',fileAbs{1});end
    end
    filePhase=strcat(fil,suff{n},'Ph',outExt);
    fileMatl=strcat(fil,suff{n},'.mat');
    fileMatlRaw=strcat(fil,suff{n},'_raw','.mat');
                
    x{n}=single(nii.img);
    if gpu;x{n}=gpuArray(x{n});end
    if exist(filePhase,'file')
        if headonly;nii.hdr=load_untouch_header_only(filePhase);nii.img=[];else nii=load_untouch_nii(filePhase);end
        p=nii.img;
        if gpu;p=gpuArray(p);end
        x{n}=x{n}.*exp(1i*p);
    end
    if exist(fileMatl,'file')
        load(fileMatl,'dat');%Dat is prescribed name for metadata
        datima{n}=dat;
    end        
    if exist(fileMatlRaw,'file')
        load(fileMatlRaw,'Par');%Par is prescribed name for rawdata
        datima{n}.Raw=Par;
    end
    hdrU{n}=nii.hdr;
    MS{n}=hdrU{n}.dime.pixdim(2:4);
    
    MT{n}=eye(4);MT{n}(1,:)=hdrU{n}.hist.srow_x;MT{n}(2,:)=hdrU{n}.hist.srow_y;MT{n}(3,:)=hdrU{n}.hist.srow_z;
    if (all(hdrU{n}.hist.srow_x==0) && all(hdrU{n}.hist.srow_y==0) && all(hdrU{n}.hist.srow_z==0)) || hdrU{n}.hist.qform_code==1%%SECOND CONDITION IS HIGHLY EXPERIMENTAL!!
        T=eye(4,4);
        T(1:3,4)=[hdrU{n}.hist.qoffset_x hdrU{n}.hist.qoffset_y hdrU{n}.hist.qoffset_z];
        q=[hdrU{n}.hist.quatern_b hdrU{n}.hist.quatern_c hdrU{n}.hist.quatern_d];
        nq=norm(q);if nq>1;q=q/nq;end;nq=min(nq,1);%To avoid complex numbers due to rounding
        q=[sqrt(1-nq^2) q];
        %if abs(hdrU{n}.hist.quatern_c)>abs(hdrU{n}.hist.quatern_d);q(1)=-q(1);end
        R=eye(4);R(1:3,1:3)=quat2rotm(q);
        S=ones(1,4);S(1:3)=MS{n};S=diag(S);
        MT{n}=T*R*S;
        
        qfac=hdrU{n}.dime.pixdim(1);
        if qfac==0;qfac=1;end        
        if qfac==-1
            D=eye(4);D(3,3)=-1;
            x{n}=flip(x{n},3);                        
            T(:,4)=MT{n}*D*[0;0;size(x{n},3)-1;1];           
            MT{n}=T*R*S;
        end                              
    end
end
if N==0;x=x{1};MS=MS{1};MT=MT{1};datima=datima{1};hdrU=hdrU{1};end