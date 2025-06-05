
function R = Rx (thetax)

%RX defines the transformation matrix for a given rotation around the x-axis.
%   [R]=RX(THETAX)
%   * THETAX is the rotation around the x-axis.
%   ** R is the transformation matrix.
%
%   Yannick Brackenier 2022-01-30

R = eye(4,'like', thetax ); 
R(1:3,1:3) = [[ 1     , 0             , 0           ]; ...      
              [ 0     , cos(thetax)   ,-sin(thetax) ]; ...
              [ 0     , sin(thetax)   ,cos(thetax)  ]];
end