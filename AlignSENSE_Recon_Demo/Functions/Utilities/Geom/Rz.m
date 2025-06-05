
function R = Rz (thetaz)

%RZ defines the transformation matrix for a given rotation around the z-axis.
%   [R]=RZ(THETAZ)
%   * THETAZ is the rotation around the z-axis.
%   ** R is the transformation matrix.
%
%   Yannick Brackenier 2022-01-30

R = eye(4,'like', thetaz ); 
R(1:3,1:3) = [[ cos(thetaz)   , -sin(thetaz)  ,0 ]; ...         
              [ sin(thetaz)   , cos(thetaz)   ,0 ]; ...
              [ 0             , 0             ,1 ]];
end
