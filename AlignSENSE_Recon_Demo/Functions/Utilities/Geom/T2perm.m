

function [perm, fl] = T2perm(MT)

%T2PERM finds a flip and permution order that (when applied) approximates a given transformation matrix.
%   [PERM,FL]=T2PERM(MT)
%   * MT the transformation matrix.
%   ** PERM the permutation order.
%   ** FL the flip indices for every dimension.
%
%   Yannick Brackenier 2022-01-30


MS=sqrt(sum(MT(1:3,1:3).^2,1));
MT = MT(1:3,1:3);%Extract only 3 dimensions, but can be extended to N dimensions

R = MT/diag(MS);% Same as MT*inv(diag(MS));
perm = 1:size(R,1);
fl = zeros(1,size(R,1));

[~,idx] = max(abs(R),[],1);
[~,idxMin] = min((R),[],1);%Since Tflip multiplies right, the columns are multiplied with -1, so need to check in that dimension

perm(idx) = perm;
fl(idx==idxMin) = 1;

%first flip and then permute
%So T = Tpermute * Tflip where Tflip and Tpermute are defined around FIXED axes
%So when in flipPermute you want to express an array in NIFTI frame, you
%have to apply forward!


%If you have an array and a corresponding MT, and you want to know which
%logical array is the Foot Head direction
%[permTemp,flTemp] = T2perm(MT);dimFH = permTemp(3);
%To know which direction: flTemp(dimIS)==1 means that inverted : logical
%axis corresponds to HF