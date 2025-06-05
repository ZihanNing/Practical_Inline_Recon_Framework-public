function Pg=phaseMapGradient(PP,dims,N)

%PHASEMAPGRADIENT   Builds a phase map gradient from linear phase 
%parameters
%   PG=PHASEMAPGRADIENT(PP,DIMS,N)
%   * PP are a set of phase parameters
%   * DIMS are the dimensions of the linear terms
%   * N are the sizes of the array on which this map is to be applied
%   * PG is the generated map as a cell over the different parameters
%

gpu=isa(PP,'gpuArray');

LD=length(dims);
assert(size(PP,dims(1))-1==LD,'Number of linear terms (%d) not matched with number of dimension on which to apply the map (%d)',size(PP,dims(1))-1,LD);
N(end+1:max(dims))=1;
NG=ones(1,length(N));NGC=NG;NG(dims)=N(dims);
rGrid=generateGrid(NG,gpu,NG,NGC);
P=dynInd(PP,1,dims(1));
for n=1:LD;P=bsxfun(@plus,P,bsxfun(@times,rGrid{dims(n)},dynInd(PP,n+1,dims(1))));end
P=1i*exp(1i*P);
Pg=cell(1,LD+1);
Pg{1}=P;
for n=1:LD;Pg{n+1}=bsxfun(@times,P,rGrid{dims(n)});end
