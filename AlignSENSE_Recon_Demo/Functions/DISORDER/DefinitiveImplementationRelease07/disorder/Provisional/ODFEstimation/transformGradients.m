function gr=transformGradients(gr,T,di)

%TRANSFORMGRADIENTS   Transforms the DWI gradients using a set of motion
%estimates
%   GR=TRANSFORMGRADIENTS(GR,T,{DI})
%   * GR are the gradients to be transformed
%   * T are the motion estimates
%   * {DI} is the transform direction, generally forward (1)
%   * GR are the transformed gradients
%

if ~exist('di','var') || isempty(di);di=1;end

gpu=isa(gr,'gpuArray');

tr=[1 3 2;
    2 1 3
    3 2 1];
ND=numDims(T);NT=size(T);

if di;mV=1:3;else mV=3:-1:1;end
for m=1:mV
    R=single(zeros([3 3 NT(3:ND-1)]));
    if gpu;R=gpuArray(R);end
    th=dynInd(T,m+3,ND);if ~di;th=-th;end
    R=dynInd(R,{tr(1:2,m),tr(1:2,m)},1:2,vertcat(horzcat(cos(th),-sin(th)),horzcat(sin(th),cos(th))));
    R=dynInd(R,[tr(3,m) tr(3,m)],1:2,1);
    gr=matfun(@mtimes,R,gr);
end
gr=reshape(gr,[3 1 NT(3:ND-1)]);