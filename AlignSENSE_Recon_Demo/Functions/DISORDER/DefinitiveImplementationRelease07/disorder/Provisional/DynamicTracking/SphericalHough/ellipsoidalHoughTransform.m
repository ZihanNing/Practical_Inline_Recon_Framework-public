function [x,sphcen,sphrad]=ellipsoidalHoughTransform(x,voxsiz,rRange,outA,grTh,Mask)

%ELLIPSOIDALHOUGHTRANSFORM detects spherical structures from 3D images. 
%Determines object centers and radii and outputs image mask of centers and 
%spheres. The Hough transform is based on the gradient field of the image.
%This code inherits from: http://www.mathworks.com/matlabcentral/fileexchange/48219
%which in turn inherits from:
%   2007_03_07 original Circular Hough Transform by Yao Peng
%   2010_08_25 Spherical Hough Transform by Brian Hulette
%   2014_10_20 simpification, filtration, image mapping, classifiers, and 
%              example by Luke Xie
%Actually it is a modification of previous spherical detection using the
%radius of curvature of the image gradients
%   X=ELLIPSOIDALHOUGHTRANSFORM(X,VOXSIZ,RRANGE,{OUTA},{GRTH},{MEDFILTSIZE},{MASK})
%   * X is the data on which to compute the transform
%   * VOXSIZ is the voxel size
%   * RRANGE is the range of valid radius in mm
%   * {OUTA} is a flag encompassing different options for function returns,
%   see the code for the details
%   * {GRTH} is the thresholding of the gradient magnitude to be performed
%   before the voting process of the transform, i.e., pixels with gradient
%   magnitudes smaller than 'grdthres' are not consider in the computation.
%   Defaults to 0.2, has to be between 0 and 1.
%   * {MASK} is a mask to constrain the detection to foreground voxels
%   * X is the resulting sphere image (or other feature according to quick
%   return options)
%   * SPHCEN is the estimated center
%   * SPHRAD is the estimated radius
%

%DEFAULT VALUES
if nargin<4 || isempty(outA);outA=0;end
if nargin<5 || isempty(grTh);grTh=0.1;end%0.01;end
if nargin<6 || isempty(Mask);Mask=ones(size(x),'like',x);end

gpu=isa(x,'gpuArray');

%writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'x'},{x});


%ARGUMENT VALIDATION
assert(ndims(x)==3,'The image has to be 3D and it is %d',ndims(x));
assert(isnumeric(x),'The image has to be numeric');
assert(numel(rRange)==2,'Range of radii has to be a two-element vector and it is a %d-element vector',numel(rRange));
rRange=sort(max([0,0;rRange(1),rRange(2)]));
assert(grTh>=0 && grTh<=1,'Gradient threshold has to be in the range 0-1 and it is %.2f',grTh);
N=size(x);N(end+1:3)=1;
ND=length(N);

rRange=round(rRange/(prod(voxsiz).^(1/3)));%Radious in pixels

H=buildFilter(N(1:3),'tukeyIso',voxsiz/2,gpu,1,0);
x=abs(filtering(x,H));

%COMPUTE AND THRESHOLD THE GRADIENT
[gm,gn,si,rv]=geometricFeatures(x,grTh);

%gm>grTh*max(gm(:)) to detect edges
%si>=0.5 to detect darker to brighter edges
%rRange(1)<=rv<=rRange(2) to detect edges with admissible curvature
indGrdMask = find(gm > grTh*max(gm(:)) & si>=0.5 & rv<=rRange(2) & rv>=rRange(1));
[indGrMask{1},indGrMask{2},indGrMask{3}]=ind2sub(N,indGrdMask);

%INDICES AND SUBSCRIPTS OF ALL THE VOTINGS TO THE ACCUMULATION
%A row in matrix lin2accum contains the indexes of the votings that are
%introduced by the different voxels in the image, which depends on the 
%radius of curvature of the gradient at those voxels
lin2Acc=cell(1,3);
for n=1:3
    grx=dynInd(gn,n,ND+1);
    lin2Acc{n}=floor(bsxfun(@plus,-bsxfun(@times,grx(indGrdMask),rv(indGrdMask)),indGrMask{n}+0.5));% Compute radius offsets and add to base position
end

% Clip the votings that are out of the accumulation array
maskValid=(lin2Acc{3}>0 & lin2Acc{3}<N(3)+1 & lin2Acc{2}>0 & lin2Acc{2}<N(2)+1 & lin2Acc{1}>0 & lin2Acc{1}<N(1)+1);
maskNonValid=~maskValid;
for n=1:3;lin2Acc{n}=lin2Acc{n}.*maskValid+maskNonValid;end
clear maskNonValid

% Linear indices (of the votings) into the accumulation array
lin2acc=sub2ind(N,lin2Acc{1},lin2Acc{2},lin2Acc{3});
clear lin2Acc

% Weights of the votings, currently using the gradient magnitudes
% but in fact any scheme can be used (application dependent)
w2acc=bsxfun(@times,gm(indGrdMask),maskValid);
%w2acc=bsxfun(@times,gm(indGrdMask).*((2*(si(indGrdMask)-0.5)).^2),maskValid);
clear maskValid
% Build the accumulation array using Matlab function 'accumarray'

accum=accumarray(lin2acc(:) ,w2acc(:));
clear w2acc
accum(end+1:prod(N))=0;
accum=reshape(accum,N);
accum=accum.*x;%To prevent clusters in too dark areas
accum=accum.*Mask;%Only voxels in the foreground are considered

if outA==3;sphcen=[];sphrad=[];x=accum;return;end

%writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'acc'},{accum});

%%LOCATE CENTERS
%H=buildFilter(N,'tukeyIso',voxsiz/16,gpu,1);
H=buildFilter(N,'tukeyIso',voxsiz/32,gpu,1);%RECENT CHANGE
accum=abs(filtering(accum,H));

%writeNII('/home/lcg13/Work/DataDefinitiveImplementationDebug06/x',{'accF'},{accum});

if outA==1;sphcen=[];sphrad=[];x=accum;return;end

[maxAccum,indMaccum]=max(accum(:));
[sphcen(1),sphcen(2),sphcen(3)]=ind2sub(N,gather(indMaccum));
fprintf('Estimated center:%s\n',sprintf(' %d',sphcen));

if rRange(1)==rRange(2) || outA==2
    %BUILD RESULTING MASK IMAGE
    sphrad=rRange(1);
    xGrid=generateGrid(N,gpu,N,[0 0 0],-sphcen);
    rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));
    x=single(rGrid<sphrad);
    return
end

%THE METHOD SHOULD PROBABLY BE REVIEWED AND SIMPLIFIED FROM HERE ON, BUT IT
%IS ALREADY QUITE ROBUST. A POTENTIAL IMPROVEMENT COULD COME FROM COMPARING
%THE DISTANCES TO THE CENTER WITH THE RADIOUS OF CURVATURE, AND DECIDE ON THE
%MOST APPROPRIATE RADIOUS FROM THEM

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
for n=1:3
    grx=dynInd(gn,n,ND+1).*gm;
    f4SgnCv=f4SgnCv+bsxfun(@times,dynInd(grx,SgnCvMatAux,1:3),SgnCvMat{n});
end
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

% Smooth the sgn. curve       
SgnCv=filtfilt(fltr4SgnCv,1,double(gather(SgnCv)));

% Get the maximum value in the sgn. curve
[~,sphrad]=max(SgnCv);
fprintf('Estimated radius:%s\n',sprintf(' %d',sphrad));

%BUILD RESULTING MASK IMAGE
xGrid=generateGrid(N,gpu,N,[0 0 0],-sphcen);
rGrid=sqrt(bsxfun(@plus,bsxfun(@plus,xGrid{1}.^2,xGrid{2}.^2),xGrid{3}.^2));
x=single(rGrid<sphrad);

