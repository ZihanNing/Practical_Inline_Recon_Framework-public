%%%%%%%%
%
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)


function x = unring(x, params);
    if nargin <2;params = [1 3 20];end

    x=double(x);
    x=ringRm(x, params);
    x=single(x);

