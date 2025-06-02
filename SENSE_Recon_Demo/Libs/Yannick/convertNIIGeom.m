

function [MTOu] = convertNIIGeom(MSIn, MTIn, typeIn, typeOu)

%CONVERTNIIGEOM converts the NIFTI geometry (orientation) information between s-form and q-form representation.
%               Translations are not included as they can be directly copied between s-form and q-form.
%   [MTOU,MTIN]=CONVERTNIIGEOM(MTIN,MSIN,TYPEIN,TYPEOU)
%   * MSIN the resolution associated with the input type. Only used to go from qForm to sForm as the latter includes resolution.
%   * MTIN the orienatation information of the input type.
%   * TYPEIN the input type of the orientation information.
%   * TYPEOU the output type of the orientation information.
%   ** MTOU the output orientation information.
%
%   Yannick Brackenier 2022-01-30

if isempty(MSIn)
    if strcmp( typeIn, 'sForm')
        MSIn = sqrt(sum(MTIn(1:3,1:3).^2,1));
        %Not used for sForm!
    else
        MSIn = ones([1 3]);
        warning('convertNIIGeom:: MS not provided for conversion to sForm. Returned sForm will not include affine part, so it must still be added to be a valid NIFTI s-form.')
    end
end

if strcmp( typeIn, 'qForm')
    assert(strcmp( typeOu, 'sForm'),'convertNIIGeom:: Output type for orientation information is not valid.')
    assert(size(MTIn,1)==1,'convertNIIGeom:: qForm representation is not valid. It should be a row.')

    %%% Complete quaternion information if only 3 provided.
    if size(MTIn,2)==3; MTIn = cat(2, sqrt( 1 -  ( MTIn(1).^2 + MTIn(2).^2 + MTIn(3).^2)),MTIn);end
    
    %%% Compute rotation matrix
    MTOu = SpinCalc('QtoDCM',[MTIn(2:4) MTIn(1)],.1,1)'; %Transpose because of DCM convention in SpinCalc
    
    %just a test
    MTOuTemp =  quaternionToR(MTIn);%Takes as input [1 3] as well as [1 4] sized arrays
    %assert( all(all( abs((MTOuTemp-MTOu)./(MTOu+eps(class(MTOu))))<eps('single')  )),'test');
    assert( all(all( abs((MTOuTemp-MTOu)./(MTOu+eps(class(MTOu))))<1e-3  )),'test');
    
    %%% Include resolution
    MTOu = MTOu * diag(MSIn);
    
elseif strcmp( typeIn, 'sForm')
    assert(strcmp( typeOu, 'qForm'),'convertNIIGeom:: Output type for orientation iformation not valid.')
    isValidGeom(MSIn, MTIn,2);
    
    R = MTIn(1:3,1:3)/diag(MSIn);
    Q = SpinCalc('DCMtoQ',R',.1,1); %Transpose because of DCM convention in SpinCalc
    
    %%% Quaternion coefficients
    MTOu = [Q(4) Q(1:3)];%These are the quaternions to be writte to NIFTI header fields. Make sure pixdim[0] is set to 1! (see https://brainder.org/2012/09/23/the-nifti-file-format/)
    
else
    error('convertNIIGeom:: Input type is not recognised.'); 
    
end