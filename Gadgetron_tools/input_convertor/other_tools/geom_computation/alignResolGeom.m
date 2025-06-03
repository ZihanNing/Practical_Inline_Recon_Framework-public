function [newVoxelSize,APhiRec,APhiRecOrig]=alignResolGeom(newMatrixSize,oriMatrixSize,oriVoxelSize,APhiRec,APhiRecOrig,permuteHist)
%%% This is a function to convert geom based on RTarget 
%%% align the geom with new matrix size and voxel size based on RTarget
%%%
%%% INPUT:
%%% newMatrixSize: the targeted matrix size 
%%% oriMatrixSize
%%% oriVoxelSize
%%% APhiRecOrig: original geom matrix with [Lin Col Par...]
%%% AphiRec: original geom matrix with [Col Lin Par...]
%%%
%%% OUTPUT:
%%% newVoxelSize
%%% newAPhiRec: updated geom matrix with [Col Lin Par...]
%%% newAPhiRecOrig: updated geom matrix with [Lin Col Par...]
%%%
%%% by Zihan @ king's
%%% March-2025

if nargin<4 || nargin<5 || isempty(APhiRecOrig) || isempty(permuteHist);APhiRecOrig=[];permuteHist = [];end

%Change orientation information accordingly
for l=2:3 % only for PE directions
    if ~isempty(APhiRecOrig) && ~isempty(permuteHist); permutegeom = 1; end
    if permutegeom; [~,APhiRecOrig]=mapNIIGeom([],APhiRecOrig,'ipermute',permuteHist{1});end
    if oriMatrixSize(l)>newMatrixSize(l)%Extract
        vr=(centerIdx(oriMatrixSize(l))  - centerIdx(newMatrixSize(l)))+[1:newMatrixSize(l)];
        dynIndParam = {':',':',':'};dynIndParam(l) = {vr};
        [newVoxelSize , APhiRec] = mapNIIGeom(oriVoxelSize, APhiRec,'dynInd',dynIndParam,oriMatrixSize,newMatrixSize);
        if permutegeom; [~,APhiRecOrig]=mapNIIGeom([],APhiRecOrig,'dynInd',dynIndParam,oriMatrixSize,newMatrixSize);end
    else%Pad
        padDim = [0 0 0];
        padDim(l) = (centerIdx(newMatrixSize(l)) - centerIdx(oriMatrixSize(l))) ;%Only need padding before
        [newVoxelSize , APhiRec] = mapNIIGeom(oriVoxelSize , APhiRec,'padArrayND',padDim,oriMatrixSize,newMatrixSize);
        if permutegeom; [~,APhiRecOrig]=mapNIIGeom([],APhiRecOrig,'padArrayND',padDim,oriMatrixSize,newMatrixSize); end
    end
    if permutegeom; [~,APhiRecOrig]=mapNIIGeom([],APhiRecOrig,'permute',permuteHist{1}); else; APhiRecOrig = []; end
end
