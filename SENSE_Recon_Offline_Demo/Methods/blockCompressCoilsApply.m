function yo=blockCompressCoilsApply(perc,y,Ao,useFold,useInvMot,quick,rad)

% BLOCKCOMPRESSCOILSAPPLY performs channel compression based on [1] M Buehrer, 
%   KP Pruessmann, P Boesiger, S Kozerke, "Array Compression for MRI With 
%   Large Coil Arrays," Magn Reson Med, 57:1131-1139, 2007, using blocks as
%   described in [2] B Bilgic, JP Marques, LL Wald, K Setsompop, "Block 
%   Coil Compression for Virtual Body Coil without Phase Singularities," 
%   ISMRM, 2016, which are aligned as in [3] T Zhang, JM Pauly, SS 
%   Vasanawala, M Lustig, "Coil Compression for Accelerated Imaging with 
%   Cartesian Sampling," Magn Reson Med 69:571–582, 2013
%   [S,Y,C,D,U]=BLOCKCOMPRESSCOILS(S,PERC,{Y},{USEFOLD},{USEINVMOT}) 
%   * S is the original channel information
%   * PERC is the number of components to be kept if bigger or equal than 
%   one. If empty no compression is performed
%   * {Y} is the original data sampled across channels
%   * {USEFOLD} is a flag that denotes whether to use the folding 
%   structure when compressing. Defaults to 1 if {Y} is present and to {0} 
%   otherwise
%   * {USEINVMOT} is a flag that indicates to use the inverse coils for 
%   motion estimation, deprecated
%   * {RAD} is a radious for blocks, it defaults to 2 (five slices in each
%   block)
%   ** So is the compressed channel information
%   ** Yo is the compressed data
%

if nargin<4 || isempty(useFold);useFold=1;end
if nargin<5 || isempty(useInvMot);useInvMot=0;end
if nargin<6 || isempty(quick);quick=0;end
if nargin<7 || isempty(rad);rad=2;end
wi=cast(resPop(triang(2*rad+1),1,[],3),'like',y);
gpu=isa(Y,'gpuArray');gpuIn=useGPU;

NY=size(y);NY(4)=perc;yo=zeros(NY,'like',y);
wo=zeros([1 1 NY(3)],'like',y);

for sl=1+rad:NS(3)-rad
    ySl=dynInd(y,sl-rad:sl+rad,3);
    ySl=aplGPU(Ao{sl},ySl,4);%Align coils   
    yo=dynInd(yo,sl-rad:sl+rad,3,wi.*ySl+dynInd(yo,sl-rad:sl+rad,3));
    wo=dynInd(wo,sl-rad:sl+rad,3,wi+dynInd(wo,sl-rad:sl+rad,3));
end
yo=yo./wo;
