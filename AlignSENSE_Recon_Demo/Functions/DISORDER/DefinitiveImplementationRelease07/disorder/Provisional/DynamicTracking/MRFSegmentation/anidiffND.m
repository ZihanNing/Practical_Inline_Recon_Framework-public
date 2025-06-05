function x=anidiffND(x,nIt,kappa,opt,dt,voxsiz)

%ANIDIFFND performs conventional anisotropic diffusion (Perona & Malik).
%A ND neighboring system structure using the structuring element [-1 0 1]
%is considered for diffusion conduction. The method is a generalization of
%the code available at
%https://uk.mathworks.com/matlabcentral/fileexchange/14995-anisotropic-diffusion--perona---malik-?requestedDomain=true
%which is credited as
%   MATLAB implementation based on Peter Kovesi's anisodiff(.):
%   P. D. Kovesi. MATLAB and Octave Functions for Computer Vision and Image Processing.
%   School of Computer Science & Software Engineering,
%   The University of Western Australia. Available from:
%   <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
% 
% Credits:
% Daniel Simoes Lopes
% ICIST
% Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% danlopes (at) civil ist utl pt
% http://www.civil.ist.utl.pt/~danlopes
%
% May 2007 original version.
%which in turns relies on [1] P Perona and J. Malik, "Scale-Space and Edge 
%Detection Using Anisotropic Diffusion," IEEE PAMI, 12(7):629-639, Jul 1990
%and [2] G Gerig, et al, "Nonlinear Anisotropic Filtering of MRI Data,"
%IEEE TMI, 11(2):221-232, Jun 1992.
%   X=ANIDIFFND(X,NIT,DT,K,OPTION,VOXSIZ)
%   * X is the data on which to compute the anisotropic diffusion
%   * NIT is the number of iterations
%   * KAPPA is the gradient modulus threshold that controls the conduction
%   * OPT serves to use one of the conduction coefficient functions 
%   proposed by Perona & Malik:
%   1 - c(x,y,z,t) = exp(-(nablaI/kappa).^2), privileges high-contrast 
%   edges over low-contrast ones. 
%   2 - c(x,y,z,t) = 1./(1 + (nablaI/kappa).^2), privileges wide regions 
%   over smaller ones.
%   Defaults to 2
%   * DT is the integration constant. 0<1/(1+sum_n 1/(Delta_n)^2) with n the
%   number of elements of the neighboring system and Delta_n the distances
%   among elements. Defaults to the maximum value
%   * VOXSIZ indicates the distances between adjactent elements of the 
%   array in each dimension. It defaults to ones(1,ND) with ND the number
%   of dimensions of the array
%   * X is the (diffused) volume with the largest scale-space parameter


if nargin<4 || isempty(opt);opt=2;end

ND=numDims(x);
if nargin<6 || isempty(voxsiz);voxsiz=ones(1,ND);end
%CURRENTLY THE NUMBER OF DIMENSIONS HAS TO MATCH THE LENGTH OF THE VOXEL
%SPACINGS BUT IN THE FUTURE IT SHOULD BE MADE GENERIC SO AS TO DIFFUSE ONLY
%ALONG THE FIRST N DIMENSIONS AS GIVEN BY THE LENGTH OF THE VOXEL SPACINGS
assert(length(voxsiz)==ND,'Data dimensions (%d) do not match voxel spacings (%d)',ND,length(voxsiz));

Nn=3*ones(1,ND);%[-1 0 1] neighboring system
Ns=1:3^ND;%Elements within the neighboring system
neigh=ind2subV(Nn,Ns)-2;
neigh(all(neigh==0,2),:)=[];%Delete current voxel from neighboring system
NE=size(neigh,1);
dv=sqrt(sum(bsxfun(@times,abs(neigh),voxsiz.^2),2));%Spacing of the neighborhood system
dv=1./(dv.^2);

h=zeros([Nn NE],'like',x);
h=dynInd(h,2*ones(1,ND),1:ND,-1);
ih=bsxfun(@plus,2*ones(1,ND),neigh);
for n=1:NE;h=dynInd(h,[ih(n,:) n],1:ND+1,1);end
perm=1:ND+1;perm([1 ND+1])=[ND+1 1];
dv=permute(dv,perm);
if nargin<5 || isempty(dt);dt=1./sum(dv(:));end%THIS MAY BE MADE GENERAL
assert(dt<=1./sum(dv(:))+eps,'Too high integration constant. Aborted method due to potential numerical instabilities');

if isempty(kappa)
    kappa=convnfft(x,h);
    kappa=sqrt(2)*sqrt(mean(abs(kappa(:)).^2));
    fprintf('Kappa: %s\n',num2str(kappa,'%.2f'));
end

%ANISOTROPIC DIFFUSION
for t=1:nIt
    %THIS SHOULD RUN RECURSIVELY IN PRESCRIBED ROIS ON THE IMAGE TO FIT THE
    %GPU, SIMILARLY TO CDFFILT
    %nabla=repmat(x,[ones(1,ND) NE]);    
    %for n=1:NE;nabla=dynInd(nabla,n,ND+1,convn(dynInd(nabla,n,ND+1),dynInd(h,n,ND+1),'same'));end
    nabla=convnfft(x,h);
    if opt==1;nabla=nabla.*exp(-(nabla/kappa).^2);else nabla=nabla./(1+(nabla/kappa).^2);end
    x=x+dt*sum(bsxfun(@times,dv,nabla),ND+1);

    % Iteration warning.
    %fprintf('\rIteration %d\n',t);
end
