
function M = intravoxDephasing(T,c,NX,Db,Dl)

%INTRAVOXDEPHASING   Computes the signal modulation coming from intravoxel dephasing.
%   Based on "ACCELERATING ITERATIVE FIELD-COMPENSATED MR IMAGE%RECONSTRUCTION ON GPUS".
%   [M]=INTRAVOXDEPHASING(T,C,NX,DB,DL)
%   * T are the motion parameters for the different motion states.
%   * C the coefficients of the basis functions used in dephaseBasis.m
%   * NX is the size of the image
%   * DB is structure containing information for the B0 modelling with basis functions
%   * DL is structure containing information for the linear B0 model
%   ** M is the spatial magnitude modulation for the different motion states

%Db and Dl seperately inserted as need in LMsolver
NT = size(T);
NX = NX(1:3);

%%% Make field per state - units of Hz
returnField =1;
B0=zeros([NX 1 NT(5)], 'like', real(Dl.D));

if exist('Db','var'); B0 = B0 + dephaseBasis( Db.B,c,NX,Db.TE,returnField);end%In Hz
if exist('Dl','var'); B0 = B0 + dephaseRotation(dynInd(T,Dl.d,6),Dl.D,Dl.TE,returnField);end%In Hz

%%% Calculate the derivative fields per state - units of  Hz/voxel
IntraVox = zeros([NX 3 NT(5)], 'like', B0);
for dim=1:3; IntraVox = dynInd( IntraVox, dim, 4, (circshift(B0,-1,dim) - circshift(B0,+1,dim))/2 );end

% visSegment(dynInd(B0,[1,1],[4,5]),[],0)
% visSegment(dynInd(IntraVox,[1,1],[4,5]),[],0)
% visSegment(dynInd(IntraVox,[2,1],[4,5]),[],0)
% visSegment(dynInd(IntraVox,[3,1],[4,5]),[],0)
%plot_(dynInd(B0,[1,1],[4,5]) , dynInd(IntraVox,[1,1],[4,5]) , dynInd(IntraVox,[2,1],[4,5]), [], [], [], 0, 0, 12)

%%% Calculate signal loss per state 
IntraVox = Db.TE* IntraVox; %units of rotations/voxel - if 1 you have zero signal -> sincN

M = ones(NX,'like',real(Dl.D));
for dim=1:3; M = M .* sincN( dynInd(IntraVox, dim, 4));end

%plot_(dynInd(B0,[1,1],[4,5]) , dynInd(M,[1,1],[4,5]) , dynInd(M,[1,1],[4,5]), [], [0.99 1], [], 0, 0, 12)

end


%%%FUNCTIONS
function x = sincU(x)% Unnormalised sinc
    idx = (x==0);
    x(~idx) = sin(x(~idx))./x(~idx);
    x(idx)=1;
end

function x = sincN(x)% Normalised sinc
    idx = (x==0);
    x(~idx) = sin(pi*x(~idx))./(pi*x(~idx));
    x(idx)=1;
end