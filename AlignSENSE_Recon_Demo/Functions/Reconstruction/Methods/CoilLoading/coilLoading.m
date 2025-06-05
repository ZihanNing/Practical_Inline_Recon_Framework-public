
function l = coilLoading(T, vB, fun)

%COILLOADING   Computes the coil loading based on the motion parameters
%   X=DEPHASEROTATION(T,C,X)
%   * T are the motion parameters 
%   * FUN is the function handle to get loading parameters from motion parameters
%   * Y are the derivative fields
%

%%% Calculate coil loading per coil and per motion state
ndT=numDims(T);
l = fun(T);
l = multDimMea(l, 6);%Sum over motion parameters

%%% Extract relevant coils
l = dynInd(l, vB, 4);

end