function x=cellFun(x,f)

%CELLFUN   Applies a given function to the elements of a cell array
%   X=CELLFUN(X,F)
%   * X is a cell array containing the input datasets
%   * F is a function handle describing the operation to perform on those
%   datasets
%   * X is a cell array containing the output datasets
%

if iscell(x)
    for n=1:length(x);x{n}=cellFun(x{n},f);end
else
    x=f(x);
end
