
function [Tcomp, mae] = compositeTransform(T, di)

%COMPOSITETRANSFORM computes the composite transformation parameters coming
% from a sequential application of transformations. Transformation
% convention follows sincRigidTransformation.m 
%   [TCOMP, MAE]=COMPOSITETRANSFORM(T,{DI})
%   * T contains the different transformations and associated parameters in the 5th dimension
%   * {DI} determines the direction of each of the transformations. Defaults to 1 (forward).
%   ** TCOMP is the composite transform (3 translations and 3 rotations)
%   ** MAE is the MAE of the rotation matrix calculated in both was. Used for error handling.
%

NT = size(T);
NTransforms = NT(5);

if nargin <2 || isempty(di); di = ones( [1 NTransforms]);end %Remove parameter dimension

%%%Calculte composite transformation matrix
Rcomp = eye(4);
for tIdx = 1:NTransforms
    Rcomp = T_to_R (dynInd(T,tIdx,5) , di(tIdx))  * Rcomp;
end

%%% Convert to Euler angles in the specified order
tol = 0.1; ichk = 1;
euler_order = '123';% 'XYZ'; This is ZYX for rotation around fixed axis = sincRigidTransform.m convention
euler = SpinCalc(strcat('DCMtoEA', euler_order), inv(Rcomp(1:3, 1:3)), tol, ichk); % DCMtoEA converts from orthonormal rotation matrix to euler angles
% Attention: You have to take the inverse of R - see SpinCalc.m: DCM - 3x3xN multidimensional matrix which pre-multiplies a coordinate
% frame column vector to calculate its coordinates in the desired new frame. YB: This is the inverse of our definition!
    
%%% Make sure rotations are in the range [-180 180]
euler = mod(euler+180,360)-180;

%%% assign Tcomp
Tcomp = zeros(dynInd(NT, 5,2,1),'like',T);
Tcomp(1:3)= Rcomp(1:3,4);
Tcomp(4:6) = flip( convertRotation(euler,'deg','rad') ); %convert to radians / T(4:6) = z,y,x rotation

%%% test if resulting rotation matrix is close to iriginal
RcompTest = T_to_R(Tcomp,1);
mae = multDimMea(abs(Rcomp - RcompTest));
tol=1e-5;
if mae>tol;fprintf('Attention: MAE of compositeTransform conversion bigger than %d - this error might accumulate throughout algorithm\n',tol);end

end
