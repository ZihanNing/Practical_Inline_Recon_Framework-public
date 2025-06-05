function C=constructMatrixHighOrderTensorToSphericalHarmonic(order)

%This is based on the constructSetOf321Polynomials from the fanDTasia toolbox

%-fanDTasia ToolBox------------------------------------------------------------------
% This Matlab script is part of the fanDTasia ToolBox: a Matlab library for Diffusion 
% Weighted MRI (DW-MRI) Processing, Diffusion Tensor (DTI) Estimation, High-order 
% Diffusion Tensor Analysis, Tensor ODF estimation, Visualization and more.
%
% A Matlab Tutorial on DW-MRI can be found in:
% http://www.cise.ufl.edu/~abarmpou/lab/fanDTasia/tutorial.php
%
%-CITATION---------------------------------------------------------------------------
% If you use this software please cite the following work:
% A. Barmpoutis and B.C. Vemuri, "A Unified Framework for Estimating Diffusion Tensors 
% of any order with Symmetric Positive-Definite Constraints", 
% In the Proceedings of ISBI, 2010
%
%-DESCRIPTION------------------------------------------------------------------------
% This function computes the coefficients of homogenous polynomials in 3 variables of
% a given order, which correspond to squares of lower (half) order polynomials. 
%
%-USE--------------------------------------------------------------------------------
% C=constructSetOfPolynomialCoef(order);
%
% order: is the order of the computed homogeneous polynomials
% C: contains the computed list of polynomial coefficients. 
%    It is a matrix of size Mprime x (2+order)!/2(order)!
%    Mprime is defined here =321. This implementation is one way of choosing polynomials,
%    there are many different ways. A slightly different way is described in the above
%    reference (Barmpoutis et al. ISBI'10).
%
%-DISCLAIMER-------------------------------------------------------------------------
% You can use this source code for non commercial research and educational purposes 
% only without licensing fees and is provided without guarantee or warrantee expressed
% or implied. You cannot repost this file without prior written permission from the 
% authors. If you use this software please cite the following work:
% A. Barmpoutis and B.C. Vemuri, "A Unified Framework for Estimating Diffusion Tensors 
% of any order with Symmetric Positive-Definite Constraints", In Proc. of ISBI, 2010.
%
%-AUTHOR-----------------------------------------------------------------------------
% Angelos Barmpoutis, PhD
% Computer and Information Science and Engineering Department
% University of Florida, Gainesville, FL 32611, USA
% abarmpou at cise dot ufl dot edu
%------------------------------------------------------------------------------------

%fprintf(1,'Constructing matrix of 321 polynomials C...');
UnitVectors;
Mprime=321;
g=g(1:Mprime,:);
for i=0:order
    for j=0:order-i
        pop(i+1,j+1,order-i-j+1)=population(i,j,order-i-j,order);
    end
end
c=1;
for i=0:order
    for j=0:order-i
        C(:,1,c)=pop(i+1,j+1,order-i-j+1)*(g(:,1).^i).*(g(:,2).^j).*(g(:,3).^(order-i-j));
        c=c+1;
    end
end
 
THETA=atan2(g(:,1),g(:,2))';
R=sqrt(sum(g.^2,2));
PHI=acos(g(:,3)./R)';
c=1;
for l=0:2:order
    Lmn=legendre(abs(l),cos(PHI));
    corrf=sqrt((2*l+1)*factorial(l-(0:l))./(4*pi*factorial(l+(0:l))));
    Lmn=bsxfun(@times,corrf',Lmn);
    Emn=exp(1i*bsxfun(@times,(0:l)',THETA));
    for m=-l:l
        n=1+abs(m);
        if m<0
            D(:,c)=sqrt(2)*(Lmn(n,:).*real(Emn(n,:)))';
        elseif m==0
            D(:,c)=(Lmn(n,:))';
        else
            D(:,c)=sqrt(2)*((-1)^(m+1))*(Lmn(n,:).*imag(Emn(n,:)))';
        end
        c=c+1;
    end
end
C=mean(bsxfun(@times,D,C),1);
C=squeeze(C);



% a2=factorial(L-M)/factorial(L+M);
            
%             
%             Lmn(n,:)
%             
%         %l
%         %size(Lmn)
%         pause        
%     end
% end
%    
% pause
%         for m=-l:l
%             
%             
%             
% 
% if L~=0
%   Lmn=squeeze(Lmn(M+1,:,:));
% end
% 
% a1=((2*L+1)/(4*pi));
% a2=factorial(L-M)/factorial(L+M);
% C=sqrt(a1*a2);
% 
% Ymn=C*Lmn.*exp(i*M*THETA);
% 
% for k=1:length(g)
%     c=1;
%     
% 
% fprintf(1,'Done\n');