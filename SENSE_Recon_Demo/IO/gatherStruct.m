function x=gatherStruct(x)
%
%GATHERSTRUCT   Gathers all elements no matter if they are fields in a
%struct or a cell array
%   X=GATHERSTRUCT(X)
%   * X is a struct to be gathered
%   ** X is the gathered struct
%
    if isstruct(x)
        y=fieldnames(x);
        for n=1:length(y)            
            x.(y{n})=gatherStruct(x.(y{n}));
        end
    elseif iscell(x)
        x=gatherCell(x);
    elseif isnumeric(x)
        %if isa(x,'gpuArray')
            x=gather(x);
        %end
    end
end

function x=gatherCell(x)
%
%GATHERSTRUCT   Gathers all elements of a cell
%   X=GATHERCELL(X)
%   * X is a cell to be gathered
%   ** X is the gathered cell
%
    for n=1:length(x)
        if iscell(x{n})
            x{n}=gatherCell(x{n});
        elseif isnumeric(x{n})
            %if isa(x,'gpuArray')
                x{n}=gather(x{n});
            %end
        end
    end
end