function B = fraccircshift(A,shiftsize)
%fraccircshift expands circshift to fractional shifts values, using linear
%interpolation. In contrast to other approaches to non-integer shifts 
%of matrices on the base of fft2 or interp2, the number of dimensions of A 
%is not limited, and the syntax of circshift applies. For integer elements 
%of shiftsize, fracircshift and circshift give the same results. 
int = floor(shiftsize);     %integer portions of shiftsize
fra = shiftsize - int;      %fractional portions of shiftsize
dim = numel(shiftsize);
B = A;
for n = 1:numel(shiftsize)  %The dimensions are treated one after another.
    intn = int(n);
    fran = fra(n);
    shift1 = zeros(dim,1);
    shift1(n) = intn;
    shift2 = zeros(dim,1);
    shift2(n) = intn+1;
    %Linear intepolation:
    B = (1-fran)*circshift(B,shift1) + fran*circshift(B,shift2);
end
