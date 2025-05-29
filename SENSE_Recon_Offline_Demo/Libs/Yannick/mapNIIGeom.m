
function [MSOu,MTOu] = mapNIIGeom( MSIn, MTIn, operation, operationParam, NIn, NOu)

%MAPNIIGEOM   Maps (changes) the NIFTI header (orientation and resolution) to a corresponding operation (resampling, flipping, permuting, padding, dynInd).
%   [MTOU,MTIN]=MAPNIIGEOM(MTIN,MSIN,OPERATION,{OPERATIOPARAM},NIN,NOU)
%   * MTIN the orienatation information of the input array.
%   * MSIN the resolution of the input array.
%   * OPERATION a string with the operation used between input and output arrays. Current implementation has following operations implemented:
%               resampling,flip,permute,ipermute,padArrayND,dynInd,translate,rotate.
%   * {OPERATIONPARAM} an optional additional parameter that defines the operation. Look at the code below to understand the meaning for different operation parameters.
%   * NIN the size of the original array.
%   * {NOU} the size of the array after applying the operation.
%   ** MTOU the output orientation information.
%   ** MSOU the output resolution.
%
%   Yannick Brackenier 2022-01-30

if isempty(MSIn); MSIn=sqrt(sum(MTIn(1:3,1:3).^2,1));end
if nargin < 4 || isempty(operation); operation=[];end
if nargin < 4 || isempty(operationParam); operationParam=[];end
if (nargin < 5 || isempty(NIn)) && any( cat( 2, strcmp(operation,'permute'), strcmp(operation,'ipermute'),  strcmp(operation,'padArrayND'),strcmp(operation,'translate') ));NIn=[1 1 1];end
if nargin < 6 || isempty(NOu); NOu=NIn;end
NIn = NIn(1:3);
NOu = NOu(1:3);

%%% Test if NIFTI header is valid
isValidGeom(MSIn,MTIn,2);

%%% RESAMPLING: resampling an array to a different resolution
if strcmp(operation, 'resampling')

    MSOu = MSIn.* NIn./NOu;
    MTOu = eye(4);
    
    %%% Change offset because of changed voxel size (NIFTI convention that world frame is for the voxel centre)
    voxelShift = ((MSOu-MSIn)/2 ) ./ MSIn;
    voxelShift = zerosL(voxelShift);%Apparantly this shift is not needed when using the DFT resampling
    MTOu(1:3,4) = MTIn(1:3,4) + MTIn(1:3,1:3)*voxelShift';
        
    %%% Change orientation
    MTOu(1:3,1:3) = MTIn(1:3,1:3)/(diag(MSIn));%Take out voxel size
    MTOu(1:3,1:3) = MTOu(1:3,1:3)*diag(MSOu); %Insert new voxel size

%%% IMRESIZE3: Matlab command to resample a 3D-array (without the DFT) - causes slight shift difference with the RESAMPLING.m method
elseif strcmp(operation, 'imresize3')

    MSOu = MSIn.* NIn./NOu;
    MTOu = eye(4);
    
    %%% Change offset because of changed voxel size (NIFTI convention that world frame is for the voxel centre)
    voxelShift = ((MSOu-MSIn)/2 ) ./ MSIn;%Don't set the shift to zero
    MTOu(1:3,4) = MTIn(1:3,4) + MTIn(1:3,1:3)*voxelShift';
        
    %%% Change orientation
    MTOu(1:3,1:3) = MTIn(1:3,1:3)/(diag(MSIn));%Take out voxel size
    MTOu(1:3,1:3) = MTOu(1:3,1:3)*diag(MSOu); %Insert new voxel size
    
%%% FLIP: flipping an array across a specified set of dimensions.
elseif strcmp(operation, 'flip')
    assert(isequal(NIn, NOu), 'mapNIIGeom:: flip operation should not change the array size.');
    
    dimFlip = operationParam;%The dimension along which flips are performed. Size [1 numDim].
    
    MSOu = MSIn; MTOu = MTIn;
    for i=1:length(dimFlip); MTOu(:,dimFlip(i)) = -MTOu(:,dimFlip(i)) ;end%s_form modification (sign difference)
    for i=1:length(dimFlip); MTOu(:,4) = MTOu(:,4) - (NOu( dimFlip(i))-1) .* MTOu(:,dimFlip(i)) ;end %have to add a different offset since voxel indexing does not start at the center FOV.

%%% PERMUTE
elseif strcmp(operation, 'permute')
    perm = operationParam;%Permutation order
    
    MTOu = MTIn;
    MTOu(:,1:3)=MTOu(:,perm(1:3));
    MSOu =  MSIn(perm(1:3));%Change voxel spacing
    
%%% IPERMUTE
elseif strcmp(operation, 'ipermute')
    perm = operationParam;%Permutation order
    
    MTOu = MTIn;
    MTOu(:,perm(1:3)) = MTOu(:,1:3);
    MSOu(perm(1:3)) =  MSIn;%Change voxel spacing

%%% DYNIND: Dynamic indexing as implemented in dynInd.m    
elseif strcmp(operation, 'dynInd')
    idx = operationParam(1:3); %Must be cell array with indices - at least provided for the first 3 dimensions
    
    %Check if indices are consistent with array sizes provided
    numsExtracted = cellfun(@(x) length(x) , idx);
    dimsExtracted = cellfun(@(x) ~strcmp(x,':') , idx);
    assert(isequal( numsExtracted(dimsExtracted) , NOu(dimsExtracted)) , 'mapNIIGeom:: dynInd not correctly defined.');
    for i=find(~dimsExtracted); idx{i} = 1:NIn(i);end%replace ':' with indices
    
    % Check that indices don't have interleaved samples
    diff = cellfun(@(x)gradient(x),idx, 'UniformOutput' ,false);
    diffConsistent = cellfun( @(x) isequal(x, ones(size(x))), diff);
    assert(isequal( diffConsistent(dimsExtracted) , ones( [1,sum(dimsExtracted)]) ), 'mapNIIGeom:: dynInd not supported when idx with increments > 1. Can be extended in this function if wanted.');
    
    orig = [idx{1}(1); idx{2}(1); idx{3}(1)] - 1;
    MTOu = MTIn;
    MTOu(1:3,4) = MTIn(1:3,4) + MTIn(1:3,1:3)*orig;%Change origin
    MSOu = MSIn;

%%% PADARRAYND: Array padding. Only interested in the pading close to the origin so this function also works for padding that is performed at only one side of the array.     
elseif strcmp(operation, 'padArrayND')
    padDim=operationParam;if isempty(padDim); padDim = [0 0 0];end%Number of padded voxels along first 3 dimensions
    
    MSOu = MSIn;
    
    orig = -padDim';%Not -1 since zero pad should not change anything
    MTOu = MTIn;
    MTOu(1:3,4) = MTIn(1:3,4) + MTIn(1:3,1:3)*orig;%Change origin    
    
%%% TRANSLATION: When array is translated.    
elseif strcmp(operation, 'translate')
    tranPar = operationParam;%translations in all 3 dimensions in units of voxels
    
    MSOu = MSIn;
    Rtran = eye(4); Rtran(1:3,4) = reshape(tranPar,[3,1]);
    MTOu = MTIn/Rtran;%Same as MTIn*inv(Rtran)
    %Equivalent to the below
    %MTOu = MTIn;
    %MTOu(1:3,4) = MTOu(1:3,4) - MTOu(1:3,1:3)*reshape(tranPar,[3,1]);

%%% ROTATION: When array is rotated around the centre of FOV.        
elseif strcmp(operation, 'rotate')
    assert(isequal(NIn, NOu), 'mapNIIGeom:: rotate operation should not change array size.');
    
    R = operationParam(1:3,1:3);%rotation matrix for the volume defined around the center of FOV
    
    %%% Make rotation matrix so that correctly defined
    R = cat(1,R,[0 0 0]);R=cat(2,R,[0;0;0;1]);
    Rcentre = eye(4); Rcentre(1:3,4) = + ( NIn./2);%In voxel units since in moving frame
    Rrot = Rcentre*R/Rcentre;
    
    MTOu = MTIn/Rrot;%Same as MTIn*inv(Rrot) - inverse since MT goes from world --> object. By rotating the object in the array, it is like backwards rotating the frame so that the object is still in the same physical position.

    %%% Deduce resolution
    MSOu = sqrt(sum(MTOu(1:3,1:3).^2,1));
    %NOTE: this strategy would work for affine matrices as well. 
    
elseif isempty(operation)
    MSOu = MSIn;
    MTOu = MTIn;
    warning('mapNIIGeom:: Operation provided is empty. Input is returned.');
    
else
    error('mapNIIGeom:: Operation "%s" not implemented.',operation)
end
