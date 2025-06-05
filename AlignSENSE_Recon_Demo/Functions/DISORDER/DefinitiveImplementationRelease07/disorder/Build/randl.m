function x=randl(N)

%RANDL generates random numbers with Laplace distribution. The code
%follows:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Random Numbers with Laplace distribution %  
%             with MATLAB Implementation                 %
%                                                        %
% Author: M.Sc. Eng. Hristo Zhivomirov          05/01/15 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.mathworks.com/matlabcentral/fileexchange/53397-generation-of-random-numbers-with-laplace-distribution-with-matlab-implementation
%The generated array has mu=0 and sigma=1. See the code below
%   X=RANDL(N)
%   * N is the size of the array to be generated
%   * X is the generated array
%

Nw=prod(N);
u1=rand(1,Nw);
u2=rand(1,Nw);
x=log(u1./u2);
x=x-mean(x);
x=x./std(x);
x=reshape(x,N);

% function x = randl(m, n)
% % function: x  = randl(m, n)
% % m - number of matrix rows
% % n - number of matrix columns
% % x - matrix with Laplacian distributed numbers 
% %     with mu = 0 and sigma = 1 (columnwise)
% % generation of two i.i.d. sequences
% u1 = rand(m, n);
% u2 = rand(m, n);
% % generation of a matrix with Laplacian
% % distributed numbers (columwise)
% x = log(u1./u2);
% x = bsxfun(@minus, x, mean(x));
% x = bsxfun(@rdivide, x, std(x));
% end