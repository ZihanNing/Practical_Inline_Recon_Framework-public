

function x = fillCell(x, value)

%FILLCELL takes an input cell structure and fills it with a specified value
%
%   [X]=ONESL(X,{VALUE})
%   * X is the input cell array.
%   * {VALUE} the value to fill all cell elements with.
%   ** X the output cell with the filled entries.
%

if nargin < 2 || ( ~ischar(value) && isempty(value) );value=0; end %Default filling value is zero (numeric)

if ~iscell(x);x=value;return;end

for n=1:length(x)
    if iscell(x{n}); x{n} = fillCell(x{n}, value);%Recursive call to fill the 
    else
       x{n}=value; %Fill
    end
end
