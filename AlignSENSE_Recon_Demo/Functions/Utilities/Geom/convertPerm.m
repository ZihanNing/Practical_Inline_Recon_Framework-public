
function [y1,y2] = convertPerm(x)

%CONVERTPERM converts a permutation order to and back from a corresponding transformation matrix.
%   [Y1,Y2]=CONVERTPERM(X)
%   * X is the permutation order or a transformation matrix.
%   ** Y1/Y2 are the converted permutation parameters. Look at code below to understand their meaning depending on the applied conversion.
%
%   Yannick Brackenier 2022-01-30

if size(x,1)>1; di=0;else;di=1;end %di determines to go from permutation order to a transformation matrix (1), or the other way around (0).
N = size(x,2);

if di==1
    idx = 1:N;
    perm(length(perm)+1:4)=idx(length(perm)+1:4);

    %%% Forward meaning (analogy with permute.m)
    y1 = diag(ones(1,4));
    y1 (:,perm) = y1;

    %%% Inverse (analogy with ipermute.m)
    y2 = diag(ones(1,4));
    y2= y2(:,perm);

else
    Nx = size(x);
    if Nx(1)>3 || NT(2)>3; x = x(1:3,1:3);end

    [U,~,V]=svd(x);
    x = U*V'; %Taking out affine component

    perm = 1:size(x,1);
    fl = zeros(1,size(x,1));

    [~,idx] = max(abs(x),[],1);
    [~,idxMin] = min((x),[],1);

    perm(idx) = perm;
    fl(idx==idxMin) = 1;
    y1=perm;y2=fl;
    
    %first flip and then permute
    %So T = Tpermute * Tflip where Tflip and Tpermute are defined around FIXED axes
    
end
