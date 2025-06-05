function [that,amse,Rhat]=weightedShrinkage(lambda,U,V,gamma,mum2,mun2,ca)

%WEIGHTEDSHRINKAGE  applies a optimal shrinkage based on a weighted
%Frobenius loss, as described in [1] W Leeb, "Matrix denoising for 
%weighted loss functions and heterogeneous signals", arXiv, 1902.0947v1, 
%2019
%   [THAT,AMSE]=WEIGHTEDSHRINKAGE(LAMBDA,U,V,WM,WN,M,N)
%   * LAMBDA are the observed singular values
%   * U are the observed left singular (column) vectors after weighting
%   * V are the observed right singular (column) vectors after weighting
%   * GAMMA is the aspect ratio of the original problem
%   * MUM2 is the normalized squared trace of the left hand side operator
%   * MUN2 is the normalized squared trace of the right hand side operator
%   * CA is the case. 0-> unweighted, 1-> asymptotically orthogonal, 2->
%   general
%   ** AMSE is the estimated asymptotic mean square error
%   ** THAT are the shrunk singular values
%   ** RHAT is an estimate of the rank
%

%RANK ESTIMATION
Rhat=sum(lambda>(1+sqrt(gamma)));
that=lambda;that(Rhat+1:end)=0;
lambda=lambda(1:Rhat);

%POPULATION EIGENVALUE ESTIMATION
auxa=lambda.^2-1-gamma;
t2=(auxa+sqrt(max(auxa.^2-4*gamma,0)))/2;
t=sqrt(t2);

%SQUARED COSINES BETWEEN POPULATION AND SAMPLE SINGULAR VECTORS
auxa=(1-gamma./t2.^2);
cm2=auxa./(1+gamma./t2);%Rows
cn2=auxa./(1+1./t2);%Columns
DinvC=sqrt(cm2).*sqrt(cn2);

if ca==0%Unweighted case
    that(1:Rhat)=t.*DinvC; 
    amse=-(DinvC.^2-1).*t2;
    %amse=t2-((t2.^2-gamma).^2)./((t2+gamma).*(t2.^2+t2));
    amse=sum(max(amse,0));
else    
    %EXTRACT SIGNAL COMPONENTS
    U=U(:,1:Rhat);V=V(:,1:Rhat);

    %CORRESPONDING SQUARED SINES
    sm2=1-cm2;
    sn2=1-cn2;

    %SAMPLE SINGULAR VECTORS UNDER WEIGHTING
    Ua2=abs(dot(U,U,1));Ua=sqrt(Ua2);
    U=bsxfun(@rdivide,U,Ua);
    Va2=abs(dot(V,V,1));Va=sqrt(Va2);
    V=bsxfun(@rdivide,V,Va);    

    %ENERGY OF WEIGHTED SINGULAR VECTORS
    auxbm=sm2*mum2;
    alpham2=(Ua2-auxbm)./cm2;
    auxcm=cm2.*alpham2;
    auxdm=auxcm+auxbm;
    auxem=sqrt(auxcm./auxdm);
    auxbn=sn2*mun2;
    alphan2=(Va2-auxbn)./cn2;  
    auxcn=cn2.*alphan2;
    auxdn=auxcn+auxbn;
    auxen=sqrt(auxcn./auxdn);
    
    if ca==1%Asymptotically orthogonal with respect to the weights, either one sided or random signal vectors
% % %THIS AND THE RECOMPUTATION OF RHAT SERVES TO SOLVE NUMERICAL PROBLEMS WITH
% % %SINGLETONS IF THE DYNAMIC RANGE OF THE DATA IS VERY BIG. IT SEEMS MORE
% % %COMPLICATED IN THE COMPUTATIONS BELOW, SO IT IS ADVISABLE TO OPERATE
% % %WITH DOUBLES
% %         %DinvC=DinvC.*max((alpham2./auxdm).*(alphan2./auxdn),0);
% %         %Rhat=sum(DinvC~=0);that(Rhat+1:end)=0;DinvC=DinvC(1:Rhat);
% %         %that(1:Rhat)=t(1:Rhat).*DinvC;
% %         %C=auxem(1:Rhat).*auxen(1:Rhat);  
% %         %amse=-(C.*DinvC-1).*t2(1:Rhat);

        DinvC=DinvC.*(alpham2./auxdm).*(alphan2./auxdn);        
        that(1:Rhat)=t.*DinvC;
        C=auxem.*auxen;
        DinvC=DinvC.*((Ua.*Va)./(sqrt(alpham2).*sqrt(alphan2)));
        amse=-(C.*DinvC-1).*t2.*alpham2.*alphan2;
        amse=sum(max(amse,0));
     else
        %INNER PRODUCTS OF SAMPLE SINGULAR VECTORS
        Dn=abs(U'*U);
        Dm=abs(V'*V);

        %INNER PRODUCTS OF POPULATION SINGULAR VECTORS
        auxa=1./auxem;
        En=Dn.*bsxfun(@times,auxa',auxa);
        En=En-diag(diag(En))+eye(Rhat,'like',En);
        auxa=1./auxen;
        Em=Dm.*bsxfun(@times,auxa',auxa);
        Em=Em-diag(diag(Em))+eye(Rhat,'like',Em);
 
        %INNER PRODUCTS BETWEEN POPULATION AND SAMPLE SINGULAR VECTORS UNDER WEIGHTING
        Cm=bsxfun(@times,Em,auxem');
        Cn=bsxfun(@times,En,auxen');
 
        %MATRICES FOR LEAST SQUARES
        C=Cm.*Cn;
        D=Dm.*Dn;
        E=Em.*En;     
        
        %LEAST SQUARES SOLUTION
        tw=t.*sqrt(alpham2).*sqrt(alphan2);
        auxa=C*tw';
        thatw=mldivide(D,auxa);
        that(1:Rhat)=thatw'./(Ua.*Va);
        %AMSE
        %M=chol(D);
        %amse=norm(M*thatw-mldivide(M',auxa))^2-tw*(C'*mldivide(D,C)-E)*tw';           
        amse=-tw*(C'*mldivide(D,C)-E)*tw';
        %Is the first term really 0??        
    end
end
