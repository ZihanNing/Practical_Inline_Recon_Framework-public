
function [RAS] = coilCentroids(S,MT)

%COILCENTROILS computes the coil centroids in the RAS frame based on simple thresholding on sensitivity profiles. 
%   [RAS]=COILCENTROILS(S,MT,)
%   * S are the sensitivity profiles (stored in the 4th dimension).
%   * MT is the s-form of S.
%   ** RAS is an array of size [3 NCoils] that contains the R(ight)-A(nterior)-S(uperior) indices in mm w.r.t to the isocenter for each coil.
%
%   Yannick Brackenier 2023-03-21

S = abs(S);
N = size(S);

%%% FIND CENTROID IN THE LOGICAL COORDINATES
IJK= ones([4 size(S,4)]);

for i=1:size(S,4)
    %Thresholding
    temp = S(:,:,:,i);
    temp (temp>.9*max(temp(:))) = 1;    %[~,idx] = max(temp(:));
    %Find peak y averaging the location of the thresholded area
    idx = find(temp(:)==1);
    [x,y,z]=ind2sub(N(1:3),idx);
    IJK(1,i) = mean(x);
    IJK(2,i) = mean(y);
    IJK(3,i) = mean(z);
end

%%% CONVERT TO RAS FRAME
RAS = MT*IJK;
RAS = RAS(1:3,:);
