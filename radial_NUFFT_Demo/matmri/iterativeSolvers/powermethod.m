function maxEig = powermethod(Ain,x0,NitMax,doAtA,xdiffThresh)
% Use the power method to find the largest eigenvalue of either the matrix
% A or the matrix A'*A. The latter case is the default. The input A can be
% a function handle.
% 
% (c) Corey Baron 2021

if nargin<1
    unitTest;
    fprintf('Unit test successful!\n')
    return;
end

if nargin<3 || isempty(NitMax)
    NitMax = 50;
end
if nargin<4 || isempty(doAtA)
    doAtA = 1;
end
if nargin<5 || isempty(xdiffThresh)
    xdiffThresh = 1e-3;
end

% Account for different ways of supplying A
if isa(Ain,'function_handle')
    A = Ain;
else
    A = @(x,transp) Asub(x,transp,Ain);
end

% Initialize
finished = 0;
nit = 0;
x = x0;
x_prev = x0;

% Iteratively find eigenvector associated with maximum eigenvalue
while ~finished
    if doAtA
        x = A(A(x,'notransp'),'transp');
    else
        x = A(x,'notransp');
    end
    x = x/norm(x(:));
    
    % Track progress
    xdiff = norm(x(:) - x_prev(:));
    x_prev = x;
    
    nit = nit+1;
    if xdiff < xdiffThresh
        finished = 1;
    end
    if nit >= NitMax
        finished = 1;
    end
end

% Use Rayleigh quotient to find eigenvalue
if doAtA
    AtAx = A(A(x,'notransp'),'transp');
else
    AtAx = A(x,'notransp');
end
maxEig = x(:)' * AtAx(:) / (x(:)'*x(:));

end

function y = Asub(x,transp,Ain)
    if strcmp(transp,'notransp')
        y = Ain*x;
    else
        y = Ain'*x;
    end
end

function unitTest
    thresh = 1e-3;
    sz = [8,5];
    Ain = randn(sz) + 1i*randn(sz);
    x0 = randn(sz(2),1);
    e = eig(Ain'*Ain);
    eTrue = max(e);
    ePowM = powermethod(Ain,x0);
    assert(abs(eTrue-ePowM)/abs(eTrue) < thresh, 'Powermethod doAtA=1 test failed.')
    %
    sz = [5,8];
    Ain = randn(sz) + 1i*randn(sz);
	Ain = Ain'*Ain;
    x0 = randn(sz(2),1);
    e = eig(Ain);
    eTrue = max(e);
    ePowM = powermethod(Ain,x0,[],0);
    assert(abs(eTrue-ePowM)/abs(eTrue) < thresh, 'Powermethod doAtA=0 test failed.')
end