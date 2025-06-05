function x=reverseB0Field(x,B0)

%REVERSEB0FIELD   Reverses the distortion and intensity weighting effect of a B0 field
%   X=REVERSEB0FIELD(X,B0) reverses the distortion using sinc interpolation
%   * X is the image before reversing the field
%   * B0 is a B0 field in "dephasing" units
%   ** X is the image after reversing the field
%

gpu=isa(x,'gpuArray');

NPE=size(x,2);
kGrid=generateGrid(NPE,gpu,NPE,ceil((NPE+1)/2));
rGrid=generateGrid(NPE,gpu);
[DFTM,DFTMH]=buildDFTM(rGrid{1}(:),kGrid{1}(:));

ND=numDims(x);ND=max(ND,4);
[NX,NB]=parUnaFun({x,B0},@size);
NX(end+1:4)=1;NB(end+1:max(4,ND))=1;
B0=repmat(B0,NX./NB);

[x,NN]=resSub(x,3:ND);
B0=resSub(B0,3:ND);

blSzN=50;
for n=1:blSzN:NN(3)
    vN=n:min(n+blSzN-1,NN(3));
    xn=dynInd(x,vN,3);Bn=dynInd(B0,vN,3);
    DFTMB0=bsxfun(@times,permute(DFTM,[3 2 4 1]),exp(bsxfun(@times,-2*pi*1i*permute(kGrid{1},[2 3 4 1]),Bn)));
    xn=sum(bsxfun(@times,xn,DFTMB0),2);
    xn=permute(xn,[1 4 3 2]);
    xn=aplGPU(DFTMH,xn,2);
    x=dynInd(x,vN,3,xn);
end
x=reshape(x,NX);

%     
% 
% for n=1:NX(4)   
%     for m=1:NX(3)
%         xn=dynInd(x,{m,n},[3 4]);Bn=dynInd(B0,{m,n},[3 4]);
%         DFTMB0=bsxfun(@times,permute(DFTM,[3 2 4 5 6 1]),exp(bsxfun(@times,-2*pi*1i*permute(kGrid{1},[2 3 4 5 6 1]),Bn)));
%         xn=sum(bsxfun(@times,xn,DFTMB0),2);
%         xn=permute(xn,[1 6 3 4 5 2]);
%         xn=aplGPU(DFTMH,xn,2);
%         x=dynInd(x,{m,n},[3 4],xn);
%     end
% end