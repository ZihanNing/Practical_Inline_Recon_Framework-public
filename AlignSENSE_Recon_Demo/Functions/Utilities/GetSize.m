

function sizeTot = getSize(var, printSize, printThreshold)

%GETSIZE  computes the size of a variable and all the sub-variables that might be included in the variable in the form of a MATLAB structure.
%   [TOTSIZE]=GETSIZE(VAR,{PRINTSIZE},{PRINTTHRESHOLD})
%   * VAR is the variable. Can be a structure with fields nested in it.
%   * {PRINTSIZE} a flag to print the total size of the structure (1) and/or all member variables (2).
%   * {PRINTTHRESHOLD} a threshold in MB for which not to print the individual sizes of the sub-variables.
%   ** TOTSIZE the total size in bytes.
%
%   Yannick Brackenier 2022-01-30

if nargin <2 || isempty(printSize); printSize = nargout<1;end %If 2, will also get size of sub-variables
if nargin <3 || isempty(printThreshold); printThreshold =1;end %size in MB to threshold sizes. By default 1 MB.

printThreshold = printThreshold * 1e+6; %In bytes
%if printSize;fprintf('\n \n');end

%Get size
printSubVariables =  single(printSize>1);
sizeTot = getSizeGeneral(var, printSubVariables , printThreshold);

%Print
if printSize; fprintf('Total size:\n   %.2f MB / %.2f GB. \n' , sizeTot/1e6, sizeTot/1e9 );end
 
end

%%% HELPER FUNCTIONS
function sizeGeneral= getSizeGeneral(var, printSize, printThreshold)


if isstruct(var)
    sizeGeneral = 0;
    fieldNames = fieldnames(var);
    
    for i = 1:length(fieldNames)
        sizeTemp = getSizeGeneral(var.(fieldNames{i}),0,  printThreshold);%Doing it this way (calculating twice) since otherwise cannot plot the variable names.
        if (sizeTemp > printThreshold) && printSize
            fprintf('Variable size ''%s'':\n', fieldNames{i});
            sizeGeneral = sizeGeneral + getSizeGeneral(var.(fieldNames{i}),1,printThreshold );
        else
            sizeGeneral = sizeGeneral + getSizeGeneral(var.(fieldNames{i}),0 ,printThreshold);
        end
    end
else
    sizeTemp = getSizeVar(var,0);
    if (sizeTemp > printThreshold)  && printSize
        sizeGeneral = getSizeVar(var,printSize);
    else
        sizeGeneral = getSizeVar(var,printSize);
    end
end    

end

function sizeVar = getSizeVar(var,printSize)
if isa(var, 'gpuArray'); var = gather(var); end %otherwise not included 

name = getVarName(var);
size = whos(name);
sizeVar = size.bytes;

if printSize; fprintf('   %.2f MB / %.2f GB. \n' , sizeVar/1e6, sizeVar/1e9 );end
end



function out = getVarName(~)
    out = inputname(1);
end

