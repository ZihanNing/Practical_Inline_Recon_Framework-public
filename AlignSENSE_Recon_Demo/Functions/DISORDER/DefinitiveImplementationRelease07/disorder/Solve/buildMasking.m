function [C,R,discMa]=buildMasking(M,useSoftMasking,gpu)

%BUILDMASKING   Builds hard or soft masking operators
%   [C,R]=BUILDMASKING(M)
%   * M is a masking operator
%   * E is the encoding structure
%   ** XOU is the data after encoding
%

if nargin<3 || isempty(gpu);gpu=isa('M','gpuArray');end

if gpu;M=gpuArray(M);end
if useSoftMasking==2;M(:)=1;end
discMa=all(ismember(M(:),[0 1]));
NX=size(M);
if ~discMa%Soft masking
    fprintf('Using soft masking\n');
    M=abs(filtering(M,buildFilter(2*NX(1:3),'tukeyIso',1,gpu,1,1),1));
    R.Ti.la=1./(abs(M).^2+0.01);
    C=[];
else
    C.Ma=M;
    R=[];
end
