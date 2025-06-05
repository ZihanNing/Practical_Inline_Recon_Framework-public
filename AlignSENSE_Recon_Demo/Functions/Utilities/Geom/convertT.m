
function [Tou] = convertT(Tin, type, rotOrder, dimG)

%CONVERTT converts the motion parameters for the function sincRigidTranform.m. 
%   [TOU]=CONVERTT(TIN,{TYPE},{ROTORDER},{DIMG})
%   * TIN the input motion parameters
%   * {TYPE} determines the type of conversion.
%   * {ROTORDER} determines the order of Euler Angles used.
%   * {DIMG} is the dimension along which motion parameters are stored in Tin or Tou.
%   ** TOU the output motion parameters.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(type); type = 'T2R';end
if nargin < 3 || isempty(rotOrder); rotOrder = [1 2 3];end %Order of rotation around moving axis (= euler angles)
if nargin < 4 || isempty(dimG); dimG = 6;end


if strcmp(type,'T2R')
    
    ND = ndims(Tin);
    NT = size(Tin);
    assert(NT(ND-1)==1,'convertT:: Input transformation parameters expected to be one state. Currently %d.', NT(ND-1));
    
    tran = dynInd(Tin, 1:3, ND);
    rot = dynInd(Tin, 4:6, ND);%z, y and x rotation since sincRigidTransformation
    
    R = eye(4,'like',Tin);
    for i = rotOrder
        if i==1
            R = R * Rx(rot(3));
        elseif i==2
            R =  R * Ry(rot(2));
        elseif i==3
            R =  R * Rz(rot(1));
        end
    end
    R = Rt(tran) * R;
    Tou=R;
    
elseif  strcmp(type,'R2T')

    R=Tin; %Rotation matrix
    assert( abs( det(R)-1 ) <= 1e-5,'convertT :: Not a valid transformation matrix.');
    
    tol=1e-3;
    suppressWarning=1;
    eulerOrder = num2str(rotOrder);eulerOrder(eulerOrder==' ')=[];
    EA = pi/180*SpinCalc(sprintf('DCMtoEA%s', eulerOrder),R(1:3,1:3)',tol,~suppressWarning); %Transpose because of convention SpinCalc.m
    EA = mod(EA+pi, 2*pi)-pi;
    
    %on = oneStruct;
    Tou = zeros([6 1]);
    Tou(1:3) = R(1:3,4);
    Tou(4:6) = fliplr(EA)';%Flip since stored in that way in DISORDER motion parameters x-y-z should be RAS
    perm = 1:dimG; perm = cat(2,perm(2:dimG), 1);
    Tou = permute(Tou, perm);
    
else
    error('convertT:: Conversion type %s not implemented.',type)
end
    
