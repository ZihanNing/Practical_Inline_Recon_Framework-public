function y=extractSlabs(x,NZSl,ext,direc,dim)

%EXTRACTSLABS constructs a set of slabs from an image
%   Y=EXTRACTSLAB(X,NZSl,EXT,DIREC) extracts a set of slabs in z or
%   constructs the volume from the slabs
%   * X is the input image
%   * NZSl is the size of the slab
%   * EXT is a flag that determines if to extract (1) or fill (0)
%   * DIREC is a flag that determines the operator, forwards (1) or 
%   backwards (0) the slab extraction affects to
%   * DIM serves to set the dimension, it defaults to 6
%   ** Y is the output image
%

if nargin<5 || isempty(dim);dim=6;end

N=size(x);N(end+1:dim)=1;
shGr=generateGrid(NZSl,0,NZSl,ceil((NZSl+1)/2));
NZSlh=floor(NZSl/2);
shGr=shGr{1};
perm=1:14;perm([3 dim])=[dim 3];

if ext==1 && direc==1
    x=permute(x,perm);    
    y=repmat(x,[1 1 NZSl]);
    for n=1:NZSl
        Sh{dim}=-shGr(n);
        y=dynInd(y,n,3,shifting(dynInd(y,n,3),Sh));
    end
end
if ext==1 && direc==0
    x=permute(x,perm);
    y=repmat(x,[1 1 NZSl]);
    y=dynInd(y,setdiff(1:NZSl,NZSlh+1),3,0);
end
if ext==0 && direc==1
    x=dynInd(x,NZSlh+1,3);
    y=permute(x,perm);
end
if ext==0 && direc==0
    NSl=size(x,dim);
    range=zeros([NSl NZSl]);
    for s=1:NSl;range(s,:)=mod(s+shGr-1,NSl)+1;end
    
    y=zeros([N(1:2) NSl N(4:dim-1) 1],'like',x);
    for s=1:NSl;y=dynInd(y,range(s,:),3,dynInd(y,range(s,:),3)+dynInd(x,s,dim));end
end
