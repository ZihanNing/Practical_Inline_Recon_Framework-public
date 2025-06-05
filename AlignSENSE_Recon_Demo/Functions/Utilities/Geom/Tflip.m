

function [T] = Tflip(fl)

%TFLIP converts a list of dimensions to flip into the corresponding transformation matrix.
%   [T]=TFLIP(FL)
%   * FL a logical array indicating which dimension to flip.
%   ** T the corresponding transformation matrix.
%
%   Yannick Brackenier 2022-01-30

fl(end+1:3)=0;

diagElem = ones(1,4);
diagElem(fl==1)=-1;

T = diag(diagElem);