
function x=gatherStruct(x, di)

%GATHERSTRUCT is a wrapper function to recursively gather all fields in a structure.
%
%   X = GATHERSTRUCT(X,{DI})
%   * X the structure to be gatered from the GPUArray type.
%   * {DI} is the direction of the operations. 1 = gathering (Default) and 0 is for loading all fields onto the GPU.
%   ** X the gathered/gpu structure.
%      
%   Yannick Brackenier 2022-01-30

if nargin<2 || isempty(di);di=1;end

if isstruct(x)
    y=fieldnames(x);
    for n=1:length(y);x.(y{n})=gatherStruct(x.(y{n}));end
elseif iscell(x)
        x=gatherCell(x,di);
elseif isnumeric(x)
    if di; x=gather(x);else; x = gpuArray(x);end
end

end

%%% HELPER FUNCTIONS
function x=gatherCell(x,di)
    for n=1:length(x)
        if iscell(x{n})
            x{n}=gatherCell(x{n},di);
        elseif isnumeric(x{n})
            if di; x{n}=gather(x{n});else;x{n}=gpuArray(x{n});end;
        end
    end
end