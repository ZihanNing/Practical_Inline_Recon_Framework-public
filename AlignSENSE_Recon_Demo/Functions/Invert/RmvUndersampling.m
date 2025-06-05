function [data] = RmvUndersampling(data,Lin,Par,accel_R)
%%% This function is to remove the undersampled points (zero) from the
%%% kspace
%%%
%%% INPUT:
%%% data: kspace [Col, Cha, Lin, Par, Ave, Con]
%%% acq_Lin: idx of acquired points in Lin direction
%%% acq_Par: idx of acquired points in Par direction
%%% (for the above two parameters, mind to update them after you zero-padding)
%%% accel_R: accelaration factor [Lin Par]
%%% 
%%% OUTPUT:
%%% data: kspace data with undersampled points removed
%%%
%%% by Zihan @ king's
%%% March-2025

idxAcqPE = cell(1,2);
offsetAcq{1} = accel_R(1)*floor((min(Lin)-1)/accel_R(1));
offsetAcq{2} = accel_R(2)*floor((min(Par)-1)/accel_R(2));
idxAcqPE{1} = min(Lin)-offsetAcq{1}:accel_R(1):size(data,3);
idxAcqPE{2} = min(Par)-offsetAcq{2}:accel_R(2):size(data,4);

data = dynInd(data, idxAcqPE, 3:4);
fprintf('Undersampled points removed from kspace!\n');
