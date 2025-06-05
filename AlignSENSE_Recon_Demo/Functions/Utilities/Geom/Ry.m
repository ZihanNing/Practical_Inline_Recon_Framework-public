function R = Ry (thetay)

%RY defines the transformation matrix for a given rotation around the y-axis.
%   [R]=RY(THETAY)
%   * THETAY is the rotation around the y-axis.
%   ** R is the transformation matrix.
%
%   Yannick Brackenier 2022-01-30

R = eye(4,'like', thetay ); 
R(1:3,1:3) = [[ cos(thetay)   , 0         ,sin(thetay)  ]; ...   
              [ 0             , 1         ,0            ]; ...
              [ -sin(thetay)  , 0         ,cos(thetay)  ]];   
end