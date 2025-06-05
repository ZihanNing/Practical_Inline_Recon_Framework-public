function [x,r]=augmentRotation(x,only,range,cent)

%AUGMENTROTATION   Augments a given set of training images by applying 
%a random rotation
%   [X,R]=AUGMENTROTATION(X,{ONLY},{RANGE},{CENT})
%   * X is the original data
%   * {ONLY} applies the rotation only along a specific axis
%   * {RANGE} describes the range of random rotations. It defaults to 2pi
%   * {CENT} describes the center of rotation. It defaults to (N/2)+1
%   ** X is the augmented data
%   ** R is the set of applied rotations
%

if nargin<3 || isempty(range);range=2*pi;end

gpu=isa(x,'gpuArray');
N=size(x);N(end+1:5)=1;
NSamp=N(5:end);
N=N(1:3);

if nargin<4 || isempty(cent);cent=(N/2)+1;end
comp=~isreal(x);
r=zeros([ones(1,4) NSamp 6],'like',x);
NDr=numDims(r);
rot=range*(rand([prod(NSamp) 1],'like',x)-0.5);%From -pi to pi-There should be ways to apply this from -pi/4 to pi/4 by permuting axis, but that's not considered yet
if nargin>=2 && ~isempty(only)%Rotation along a specific axis    
    rot=reshape(rot,[ones(1,4) NSamp]);
    r=dynInd(r,7-only,NDr,rot);
else    
    v=randn([prod(NSamp) 3]);
    v=bsxfun(@rdivide,v,sqrt(sum(abs(v).^2,2)));%Random orientation
    rot=convertRotation(wrapTo2Pi(rot),'rad','deg');%Random rotation in degrees    
    ev=cat(2,v,rot);
    rot=SpinCalc('EVtoEA321',ev,1e-3,0);%NOTE FOR EULER ANGLES THE SINGULARITY OCCURS WHEN THE SECOND IS CLOSE TO 90, THIS IS NOT POSSIBLE HERE
    rot=convertRotation(rot,'deg','rad');
    rot=reshape(rot,[ones(1,4) NSamp 3]);
    r=dynInd(r,4:6,NDr,rot);
end
[~,kGrid,rkGrid,~,cGrid]=generateTransformGrids(N,gpu,N,cent);
et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,r,[],[],[],[],cGrid);
x=sincRigidTransform(x,et);
if ~comp;x=real(x);end
