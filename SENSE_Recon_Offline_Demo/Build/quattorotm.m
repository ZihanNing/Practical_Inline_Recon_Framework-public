function R = quattorotm( q )
%QUAT2ROTM Convert quaternion to rotation matrix
%   R = QUAT2ROTM(QOBJ) converts a quaternion object, QOBJ, into an orthonormal
%   rotation matrix, R. Each quaternion represents a 3D rotation. 
%   QOBJ is an N-element vector of quaternion objects.
%   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices.
%
%   R = QUAT2ROTM(Q) converts a unit quaternion, Q, into an orthonormal
%   rotation matrix, R. The input, Q, is an N-by-4 matrix containing N quaternions. 
%   Each quaternion represents a 3D rotation and is of the form q = [w x y z], 
%   with w as the scalar number. Each element of Q must be a real number.
%
%   Example:
%      % Convert a quaternion to rotation matrix
%      q = [0.7071 0.7071 0 0];
%      R = quat2rotm(q)
%
%      % Convert a quaternion object
%      qobj = quaternion([sqrt(2)/2 0 sqrt(2)/2 0]);
%      R = quat2rotm(qobj);
%
%   See also rotm2quat, quaternion

%   Copyright 2014-2019 The MathWorks, Inc.

%#codegen

% Validate the quaternions
%q = robotics.internal.validation.validateQuaternion(q, 'quat2rotm', 'q');

% Normalize and transpose the quaternions
%q = robotics.internal.normalizeRows(q).';
q = bsxfun(@times, q, 1./sqrt(sum(q.^2,2))).';


% Reshape the quaternions in the depth dimension
q2 = reshape(q,[4 1 size(q,2)]);

s = q2(1,1,:);
x = q2(2,1,:);
y = q2(3,1,:);
z = q2(4,1,:);

% Explicitly define concatenation dimension for codegen
tempR = cat(1, 1 - 2*(y.^2 + z.^2),   2*(x.*y - s.*z),   2*(x.*z + s.*y),...
2*(x.*y + s.*z), 1 - 2*(x.^2 + z.^2),   2*(y.*z - s.*x),...
2*(x.*z - s.*y),   2*(y.*z + s.*x), 1 - 2*(x.^2 + y.^2) );

R = reshape(tempR, [3, 3, length(s)]);
R = permute(R, [2 1 3]);

end

