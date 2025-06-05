function Bg=constructMatrixOfIntegralsV(g_, order, delta,g)
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
% This function computes the spherical integrals of the basis with all possible 
% monomials in 3 variables of a certain order. The 3 variables are given in the form 
% of 3-dimensional vectors.
%
%-USE--------------------------------------------------------------------------------
% Bg=constructMatrixOfIntegrals(g, order, delta)
%
% order: is the order of the monomials in 3 variables
% g: is a list of N 3-dimensional vectors stacked in a matrix of size Nx3
% delta: is a constant
% Bg: contains the computed basis-monomial integrals. 
%     It is a matrix of size N x (2+order)!/2(order)!
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
%fprintf(1,'Constructing matrix of Basis-Monomial Integrals...');
integr=321;
g2=g(1:integr,:);
Bg=single(zeros([size(g_,1) sum(1:order+1) size(g_,3)]));
gnorm=single(zeros([integr sum(1:order+1)]));
gpu=isa(g_,'gpuArray');
if gpu;Bg=gpuArray(Bg);gnorm=gpuArray(gnorm);end
c=1;
for l=0:order
    for m=0:order-l
        gnorm(:,c)=((g2(:,1).^l).*(g2(:,2).^m).*(g2(:,3).^(order-l-m)));
        c=c+1;
    end
end
for k=1:size(g_,1);Bg(k,:,:)=mean(bsxfun(@times,exp(-delta*sum(bsxfun(@times,g2,g_(k,:,:)),2).^2),gnorm),1);end