
function R = Rt (t)

%RT defines the transformation matrix for a given translation.
%   [R]=RT(T)
%   * T is the translation.
%   ** R is the transformation matrix.
%
%   Yannick Brackenier 2022-01-30

R = eye(4,'like',t); 
R(1:3,4) = t;

end
