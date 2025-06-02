function eul = quattoeul( q, seq )
%QUAT2EUL Convert quaternion to Euler angles
%   EUL = QUAT2EUL(QOBJ) converts a quaternion object, QOBJ, into the
%   corresponding Euler angles, EUL. Each quaternion represents 
%   a 3D rotation. QOBJ is an N-element vector of quaternion objects.
%   The output, EUL, is an N-by-3 array of Euler rotation angles with each 
%   row representing one Euler angle set. Rotation angles are in radians.
%
%   EUL = QUAT2EUL(Q) converts a unit quaternion rotation into the corresponding 
%   Euler angles. The input, Q, is an N-by-4 matrix containing N quaternions. 
%   Each quaternion represents a 3D rotation and is of the form q = [w x y z], 
%   with w as the scalar number. Each element of Q must be a real number.
%
%   EUL = QUAT2EUL(___, SEQ) converts unit quaternion into Euler angles.
%   The Euler angles are specified by the body-fixed (intrinsic) axis rotation
%   sequence, SEQ.
%
%   The default rotation sequence is 'ZYX', where the order of rotation
%   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
%
%   The following rotation sequences, SEQ, are supported: 'ZYX', 'ZYZ', and
%   'XYZ'.
%
%   Example:
%      % Calculates Euler angles for a quaternion
%      % By default, the ZYX axis order will be used.
%      q = [sqrt(2)/2 0 sqrt(2)/2 0];
%      eul = quat2eul(q)
%
%      % Calculate the Euler angles for a ZYZ rotation
%      qobj = quaternion([0.7071 0.7071 0 0]);
%      eulZYZ = quat2eul(qobj, 'ZYZ')
%
%   See also eul2quat, quaternion

%   Copyright 2014-2019 The MathWorks, Inc.

%#codegen

% Validate the quaternions
%q = robotics.internal.validation.validateQuaternion(q, 'quat2eul', 'q');

%seq = robotics.internal.validation.validateEulerSequence(varargin{:});

if nargin<2 || isempty(seq);seq='ZYX';end


% Normalize the quaternions
%q = robotics.internal.normalizeRows(q);
q = bsxfun(@times, q, 1./sqrt(sum(q.^2,2)));


qw = q(:,1);
qx = q(:,2);
qy = q(:,3);
qz = q(:,4);

% Pre-allocate output
eul = zeros(size(q,1), 3, 'like', q);

% The parsed sequence will be in all upper-case letters and validated
switch seq
    case 'ZYX'
        % Cap all inputs to asin to 1, since values >1 produce complex
        % results
        % Since the quaternion is of unit length, this should never happen,
        % but some code generation configuration seem to hit this edge case
        % under some circumstances.
        aSinInput = -2*(qx.*qz-qw.*qy);
        aSinInput(aSinInput > 1) = 1;
        aSinInput(aSinInput < -1) = -1;
        
        eul = [ atan2( 2*(qx.*qy+qw.*qz), qw.^2 + qx.^2 - qy.^2 - qz.^2 ), ...
            asin( aSinInput ), ...
            atan2( 2*(qy.*qz+qw.*qx), qw.^2 - qx.^2 - qy.^2 + qz.^2 )];
        
    case 'ZYZ'
        % Need to convert to intermediate rotation matrix here to avoid
        % singularities
        R = quattorotm(q);
        eul = rotmtoeul(R, 'ZYZ');
        
    case 'XYZ'
        % Prevent singularities as done in ZYX case
        % Alternative to rotm2eul(quat2rotm(q), 'XYZ') with fewer
        % operations
        
        % sin(y) = R13 = 2 * (qx*qz + qy*qw)
        % tan(x) = sin(x) / cos(x) = -R23 / R33
        %        = -2 * (qy*qz - qx*qw) / (1 - 2*(qx^2 + qy^2))
        %        = -2 * (qy*qz - qx*qw) / (qw^2 - qx^2 - qy^2 + qz^2)
        % tan(z) = sin(z) / cos(z) = -R12 / R11
        %        = -2 * (qx*qy - qz*qw) / (1 - 2*(qy^2 + qz^2))
        %        = -2 * (qy*qz - qx*qw) / (qw^2 + qx^2 - qy^2 - qz^2)
        
        aSinInput = 2*(qx.*qz + qy.*qw);
        aSinInput(aSinInput > 1) = 1;
        aSinInput(aSinInput < -1) = -1;
        
        eul = [ atan2( -2*(qy.*qz - qx.*qw), qw.^2 - qx.^2 - qy.^2 + qz.^2 ), ...
            asin( aSinInput ), ...
            atan2( -2*(qx.*qy - qz.*qw), qw.^2 + qx.^2 - qy.^2 - qz.^2 )];
end

% Check for complex numbers
if ~isreal(eul)
    eul = real(eul);
end

end
