function dcm2nii(fileIn,fileOu,sep)

%DCM2NII   Converts from dicom to nifti+matlab
%   DCM2NII(FILEIN,{FILEOU},{SEP})
%   * FILEIN is the input file
%   * {FILEOU} is the output file, it defaults to input file with different
%   extension
%   * {SEP} indicates whether slices are separated, it defaults to 0
%

if nargin<2 || isempty(fileOu);fileOu=removeExtension(fileIn,{'dcm'});end
if nargin<3 || isempty(sep);sep=0;end

%READ METADATA AND COMPUTE GEOMETRIC INFORMATION
MT=[];
if ~sep
    dcminf=dicominfo(fileIn);
    aux=dcminf.PerFrameFunctionalGroupsSequence.Item_1;
    %Compute spacing
    MS=[aux.PixelMeasuresSequence.Item_1.PixelSpacing' dcminf.SpacingBetweenSlices];
    %Compute transformation (THESE LINES ARE ADAPTED FROM
    %https://www.mathworks.com/matlabcentral/fileexchange/24277-gettransformmatrix)}
    ipp=aux.PlanePositionSequence.Item_1.ImagePositionPatient;
    iop=aux.PlaneOrientationSequence.Item_1.ImageOrientationPatient;
    if ~isfield(aux.PixelMeasuresSequence.Item_1,'SliceThickness');fprintf('Undefined slice thickness\n');return;end
    slTh=aux.PixelMeasuresSequence.Item_1.SliceThickness;
    computeGeometry;
    
    %READ AND ADAPT THE DATA
    x=single(dicomread(fileIn))*aux.PixelValueTransformationSequence.Item_1.RescaleSlope+aux.PixelValueTransformationSequence.Item_1.RescaleIntercept;
    N=size(x);
    aux=dcminf.PerFrameFunctionalGroupsSequence.(sprintf('Item_%d',dcminf.NumberOfFrames)).FrameContentSequence.Item_1;
    N(3:4)=[aux.InStackPositionNumber aux.TemporalPositionIndex];
    M=size(x);M(end+1:4)=1;
    if prod(M)~=prod(N);fprintf('Dimensions of data%s not compatible with dimensions in header%s: not converted\n',sprintf(' %d',M),sprintf(' %d',N));return;end
    x=reshape(x,[N(1:2) N(4) N(3)]);
    x=permute(x,[1 2 4 3]);
else
    fils=getFileNames(strcat(fileIn,'*'),1);
    N(3)=length(fils);
    for n=1:N(3);fils{n}=fullfile(fileparts(fileIn),fils{n});end

    %READ METADATA AND COMPUTE GEOMETRIC INFORMATION
    dcminf=cell(1,N(3));
    for n=1:N(3)
        dcminf{n}=dicominfo(fils{n});       
        if dcminf{n}.InstanceNumber==1
            %Compute spacing
            if isfield(dcminf{n},'SpacingBetweenSlices');MS=[dcminf{n}.PixelSpacing' dcminf{n}.SpacingBetweenSlices];%2D                
            else MS=[dcminf{n}.PixelSpacing' dcminf{n}.SliceThickness];%3D
            end
            %Compute transformation (THESE LINES ARE ADAPTED FROM
            %https://www.mathworks.com/matlabcentral/fileexchange/24277-gettransformmatrix)
            ipp=dcminf{n}.ImagePositionPatient;
            iop=dcminf{n}.ImageOrientationPatient;
            slTh=dcminf{n}.SliceThickness;
            computeGeometry;
            N(1:2)=[dcminf{n}.Width dcminf{n}.Height];
        end
    end
    
    %READ AND ADAPT THE DATA   
    NN=size(single(dicomread(fils{n})));
    if ~all(N(1:2)==NN(1:2));N([2 1])=N([1 2]);end
    x=zeros(N,'single');
    for n=1:N(3)
        if isfield(dcminf{n},'RescaleSlope');x=dynInd(x,dcminf{n}.InstanceNumber,3,single(dicomread(fils{n}))*dcminf{n}.RescaleSlope+dcminf{n}.RescaleIntercept);
        else x=dynInd(x,dcminf{n}.InstanceNumber,3,single(dicomread(fils{n})));
        end           
    end
end
x=permute(x,[2 1 3 4]);%Adapting from DICOM to NIFTI (Philips/rview combo)
dat.Dicom=dcminf;

%WRITE NIFTI
writenii(fileOu,x,[],MS,MT,dat);

function computeGeometry
    ipp
    iop
    T=[1 0 0 ipp(1); 0 1 0 ipp(2); 0 0 1 ipp(3); 0 0 0 1];%Translation
    r=iop(1:3);c=iop(4:6);s=cross(r',c');%Orientation
    R=[r(1) c(1) s(1) 0; r(2) c(2) s(2) 0; r(3) c(3) s(3) 0; 0 0 0 1];%Rotation
    S=[MS(2) 0 0 0;0 MS(1) 0 0;0 0 MS(3) 0;0 0 0 1];%Scale
    PO=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];%Permute/Flip (origin)
    PE=[-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];%Permute/Flip (this is typical for Philips/rview combo)
    MT=PE*T*R*S*PO;
    dat.Geom.PO=PO;
    dat.Geom.S=S;
    dat.Geom.R=R;
    dat.Geom.T=T;
    dat.Geom.PE=PE;
    dat.Geom.SlTh=slTh;
end

end