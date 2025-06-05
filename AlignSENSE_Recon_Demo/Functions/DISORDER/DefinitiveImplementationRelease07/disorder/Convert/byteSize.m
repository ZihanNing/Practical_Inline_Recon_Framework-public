function str=byteSize(x)

%BYTESIZE   Returns a readable string for bytes  
%code is based on https://stackoverflow.com/questions/4845561/how-to-know-the-size-of-a-variable-in-matlab
%answered Jan 31 '11 at 17:40 / MatlabSorter / 9941717
%   STR=BYTESIZE(X)
%   * X is a number of bytes
%   * STR is the returned string
%

assert(x>=0,'Number of bytes cannot be negative (%d)',x);

scale=floor(log(x)/log(1024));
switch scale
    case 0
        str=[sprintf('%.0f',x) ' b'];
    case 1
        str=[sprintf('%.2f',x/(1024)) ' kb'];
    case 2
        str=[sprintf('%.2f',x/(1024^2)) ' Mb'];
    case 3
        str=[sprintf('%.2f',x/(1024^3)) ' Gb'];
    case 4
        str=[sprintf('%.2f',x/(1024^4)) ' Tb'];
    case -inf     
        str='Not Available';% Size occasionally returned as zero
    otherwise
       str='Over a petabyte!!!';
end