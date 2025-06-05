function y=extractSlabsOld(x,NT,ext,direc,yZR)

%EXTRACTSLABS constructs a set of slabs from an image
%   Y=EXTRACTSLAB(X,NT,EXT,DIREC,YZR) extracts a set of slabs in z or
%   constructs the volume from the slabs
%   X is the input image
%   NT is the size of the slab
%   EXT is a flag that determines if to extract (1) or fill (0)
%   DIREC is a flag that determines the operator, forwards (1) or backwards
%   (0) the slab extraction affects to
%   YZR is an array with precomputed slab array sizes
%   It returns,
%   Y, the output image
%

N=size(x);
y=yZR{ext+1}{direc+1};
N(end+1:5)=1;
NTh=floor(NT/2);
NSl=N(6-3*ext);

range=zeros([NSl NT]);
for s=1:NSl
    r=s-NTh:s+NTh;
    r(r<1)=r(r<1)+NSl;
    r(r>NSl)=r(r>NSl)-NSl;
    range(s,:)=r;
end
if ext==1        
    if direc==1
        for s=1:N(3);y=dynInd(y,s,6,dynInd(x,range(s,:),3));end 
    else
        y=dynInd(y,NTh+1,3,permute(x,[1 2 6 4 5 3]));
    end
else    
    if direc==1
        y=permute(dynInd(x,NTh+1,3),[1 2 6 4 5 3]);
    else    
        for s=1:N(6);y=dynInd(y,range(s,:),3,dynInd(y,range(s,:),3)+dynInd(x,s,6));end
    end
end
