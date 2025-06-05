function [x,sphcen,sphrad]=sphericalHoughTransform(x,rRange,outA,grTh,filtLocalMax,medFiltSize)

%SPHERICALHOUGHTRANSFORM detects spherical structures from 3D images. 
%Determines object centers and radii and outputs image mask of centers and 
%spheres. The Hough transform is based on the gradient field of the image.
%This code inherits from: http://www.mathworks.com/matlabcentral/fileexchange/48219
%which in turn inherits from:
%   2007_03_07 original Circular Hough Transform by Yao Peng
%   2010_08_25 Spherical Hough Transform by Brian Hulette
%   2014_10_20 simpification, filtration, image mapping, classifiers, and 
%              example by Luke Xie
%   X=SPHERICALHOUGHTRANSFORM(X,RRANGE,{GRTH},{FILTLOCALMAX})
%   * X is the data on which to compute the transform
%   * RRANGE is the range of valid radious in pixels
%   * {GRTH} is the thresholding of the gradient magnitude to be performed
%   before the voting process of the transform, i.e., pixels with gradient
%   magnitudes smaller than 'grdthres' are not consider in the computation.
%   Defaults to 0.2, has to be between 0 and 1.
%   * {FILTLOCALMAX} is the radius of the filter used in the search of 
%   local maxima in the accumulation array. To detect spheres whose shapes 
%   are less perfect, the radius of the filter needs to be set larger.
%   Defaults to 8, minimum is 3.
%   * {MEDFILTSIZE} is the diameter of the median filter used for smoothing 
%   the accumulation array. Defaults to 7
%   * X is the resulting sphere image
%   * SPHCEN is the estimated center
%   * SPHRAD is the estimated radious
%

%DEFAULT VALUES
if nargin<3 || isempty(outA);outA=0;end
if nargin<4 || isempty(grTh);grTh=0.2;end
if nargin<5 || isempty(filtLocalMax);filtLocalMax=8;end
if nargin<6 || isempty(medFiltSize);medFiltSize=7;end


gpu=isa(x,'gpuArray');

%ARGUMENT VALIDATION
assert(ndims(x)==3,'The image has to be 3D and it is %d',ndims(x));
assert(isnumeric(x),'The image has to be numeric');
assert(numel(rRange)==2,'Range of radii has to be a two-element vector and it is a %d-element vector',numel(rRange));
rRange=sort(max([0,0;rRange(1),rRange(2)]));
assert(grTh>=0 && grTh<=1,'Gradient threshold has to be in the range 0-1 and it is %.2f',grTh);
assert(filtLocalMax>=3,'The radious of the filter to search local maxima has to be greater than 3 and it is %d',filtLocalMax);
N=size(x);

%COMPUTE AND THRESHOLD THE GRADIENT
[grx{2},grx{1},grx{3}]=gradient(x);
grMag=sqrt(grx{1}.^2+grx{2}.^2+grx{3}.^2);
indGrdMask = find(grMag > grTh*max(grMag(:)));
[indGrMask{1},indGrMask{2},indGrMask{3}]=ind2sub(N,indGrdMask);

%INDICES AND SUBSCRIPTS OF ALL THE VOTINGS TO THE ACCUMULATION
% A row in matrix 'lin2accum_aJ' contains the J indices (into the
% accumulation array) of all the votings that are introduced by a
% same pixel in the image. Similarly with matrix 'lin2accum_aI'.
if rRange(2)==rRange(1)
    dr=[-rRange(1) rRange(1)];
else
    dr =[(-rRange(2)+0.5):-rRange(1) (rRange(1)+0.5):rRange(2)];
end

%A row in matrix lin2accum contains the indices of all the votings that are
%intruced by the same pixel in the image.
lin2Acc=cell(1,3);for n=1:3;lin2Acc{n}=floor(bsxfun(@plus,(grx{n}(indGrdMask)./grMag(indGrdMask))*dr,indGrMask{n}+0.5));end% Compute radius offsets and add to base position

% Clip the votings that are out of the accumulation array
maskValid=(lin2Acc{3}>0 & lin2Acc{3}<N(3)+1 & lin2Acc{2}>0 & lin2Acc{2}<N(2)+1 & lin2Acc{1}>0 & lin2Acc{1}<N(1)+1);
maskNonValid = ~ maskValid;
for n=1:3;lin2Acc{n}=lin2Acc{n}.*maskValid+maskNonValid;end
clear maskNonValid

% Linear indices (of the votings) into the accumulation array
lin2acc=sub2ind(N,lin2Acc{1},lin2Acc{2},lin2Acc{3});
clear lin2Acc
lin2acc=lin2acc(:);

% Weights of the votings, currently using the gradient maginitudes
% but in fact any scheme can be used (application dependent)

w2acc=bsxfun(@times,grMag(indGrdMask),maskValid);
clear maskValid

% Build the accumulation array using Matlab function 'accumarray'
accum=accumarray(lin2acc ,w2acc(:));
clear w2acc
accum(end+1:prod(N))=0;
accum=reshape(accum,N);

%LOCATE CENTERS
% -- Filter for searching for local maxima
filtLocalMax_s=1.35;
filtLocalMax_r=filtLocalMax;%ceil(filtLocalMax*0.6);

%if rRange(1)==rRange(2)
%    x=accum;return;
%end

% Smooth the accumulation array using a 3D Median Filter
accum=cdfFilt(accum,'med',medFiltSize*ones(1,3));

H=buildFilter(N,'tukeyIso',ones(1,3),gpu,1);

candLM=filtering(accum,H);
candLM=abs(candLM);
if outA
x=candLM;
return
end


[maxAccum,indMaccum]=max(candLM(:));
[sphcen(1),sphcen(2),sphcen(3)]=ind2sub(N,gather(indMaccum));
fprintf('Estimated center:%s\n',sprintf(' %d',sphcen));

if rRange(1)==rRange(2);
    %BUILD RESULTING MASK IMAGE
    sphrad=rRange(1);
    xGrid=generateGrid(N,gpu,N,[0 0 0],-sphcen);
    rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));
    x=single(rGrid<sphrad);
    return
end

%% determine radii
% Parameters for the estimation of the radii of spheres
fltr4SgnCv = [2 1 1];
fltr4SgnCv = fltr4SgnCv / sum(fltr4SgnCv);

% Find sphere's radius using its signature curve

% Neighborhood region of the sphere for building the sgn. curve
circen_round = sphcen;
SCvR=single(zeros(3,2));
SgnCvMat=cell(1,3);SgnCvMatAux=SgnCvMat;
for n=1:3
    SCvR(n,1)=max(circen_round(n)-rRange(2)-1,1);
    SCvR(n,2)=min(circen_round(n)+rRange(2)+1,N(n));
    SgnCvMatAux{n}=(SCvR(n,1):SCvR(n,2));
    SgnCvMat{n}=SgnCvMatAux{n}-sphcen(n);        
end
% Sgn. curve
SgnCvMat{1}=permute(SgnCvMat{1},[2 1]);
SgnCvMat{3}=permute(SgnCvMat{3},[1 3 2]);    

SgnCvMatR=sqrt(bsxfun(@plus,bsxfun(@plus,SgnCvMat{1}.^2,SgnCvMat{2}.^2),SgnCvMat{3}.^2));
SgnCvMatRp1=round(SgnCvMatR)+1;

f4SgnCv=0;
for n=1:3;f4SgnCv=f4SgnCv+bsxfun(@times,dynInd(grx{n},SgnCvMatAux,1:3),SgnCvMat{n});end
f4SgnCv=abs(f4SgnCv)./SgnCvMatR;
SgnCv=accumarray(SgnCvMatRp1(:),f4SgnCv(:));

SgnCvCnt=accumarray(SgnCvMatRp1(:),single(ones(numel(f4SgnCv),1)));
SgnCvCnt=SgnCvCnt+(SgnCvCnt==0);
SgnCv=SgnCv./SgnCvCnt;

% Suppress the undesired entries in the sgn. curve
% -- Radii that correspond to short arcs
SgnCv=SgnCv.*(SgnCvCnt>=(pi/4*(0:(numel(SgnCvCnt)-1))'));
% -- Radii that are out of the given range
SgnCv([1:(round(rRange(1))+1) (round(rRange(2))+1):end])=0;

% Get rid of the zero radius entry in the array
SgnCv=SgnCv(2:end);

%figure;
%plot(SgnCv)

% Smooth the sgn. curve       
SgnCv=filtfilt(fltr4SgnCv,1,double(gather(SgnCv)));

%hold on
%plot(SgnCv,'--')

% Get the maximum value in the sgn. curve
[~,sphrad]=max(SgnCv);
fprintf('Estimated radious:%s\n',sprintf(' %d',sphrad));


%BUILD RESULTING MASK IMAGE
%[~,indMaxPix]=max(candrad);
xGrid=generateGrid(N,gpu,N,[0 0 0],-sphcen);
rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));
x=single(rGrid<sphrad);
