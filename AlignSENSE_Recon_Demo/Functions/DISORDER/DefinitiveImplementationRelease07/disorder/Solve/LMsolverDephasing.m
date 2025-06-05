function E=LMsolverDephasing(y,E,x,deb)

%LMSOLVERDEPHASING   Performs a Levenger-Marquardt iteration on a series of 
%   parameters of the encoding matrix to minimize a least squares 
%   backprojection of the reconstruction to the measured data
%   E=LMSOLVERDEPHASING(Y,E,X,{DEB})
%   * Y is the measured data
%   * E is the encoding structure
%   * X is the reconstructed data
%   * {DEB} indicates whether to print information about convergence
%   ** E is the modified encoding structure
%

if nargin<4 || isempty(deb);deb=2;end

%E.NMs=size(E.Tr,5);E.NSe=length(E.Fs{1});

if isfield(E,'Dc')             
    gpu=isa(x,'gpuArray');
    %GENERAL PARAMETERS AND ARRAYS
    T=E.Tr;E=rmfield(E,'Tr'); % Why removing Tr from E, can we not just let the encode(x,E) do the transformations? 
    Dc.c=E.Dc.c;Dc.D=E.Dc.D;%Coefficients we want to estimate (2 x number of atoms)
    NT=size(Dc.c);
    Sf=gather(E.Sf);
    bS=E.bS;dS=E.dS;oS=E.oS; % YB: save them as these parameters will be changed within the code to match one motion state      

    if NT(3)~=1
        dimR=3;% YB: dimension that is fully encoded
        ndR=sum(1-E.Dc.cC); %YB: number of non-converged coefficients - only in the fully encoded direction
        indRNoConv=find(~E.Dc.cC);
        y=dynInd(y,indRNoConv,3);% 3rd dimension (readout) in image space!
        E.Sf=dynInd(E.Sf,indRNoConv,3);% YB: Why extracting indices out of that dimension??? 
        Dc.c=dynInd(E.Dc.c,indRNoConv,3);
        Dc.D=dynInd(E.Dc.D,indRNoConv,3);
    else        
        dimR=[];
        ndR=1;
        Dc.c=E.Dc.c;
        Dc.D=E.Dc.D;
    end
    dimS=6:7;%Dimensions of parameters of the transform
    ndS=NT(dimS);
    ndG=prod(ndS);
    ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
    ndH=size(ha,2);% Number of combinations
    flagw=zeros([1 1 ndR]);En=flagw;EnPrev=flagw;yN=flagw;xN=flagw;%YB: why not flags per coefficients? 
    dH=single(zeros([ndH ndR]));dG=single(zeros([ndG ndR]));%Non-separable versus separable dimensions
    multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)    
    %YB: where do weights come from???
        
    nY=numel(y);
    E.bS=ceil(E.bS/ndG);
    ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;           
    %COMPUTE THE JACOBIAN AND STORE VALUES
    for a=ET.oS(1):ET.bS(1):ET.dS(1);vA=a:min(a+ET.bS(1)-1,ET.dS(1));%vA=vA(vA<=E.NMs); %YB: use ET for bS since bS in E changes within look
        indY=[];for c=1:length(vA);indY=horzcat(indY,E.nEc(vA(c))+1:E.nEc(vA(c)+1));end        
        Tr=dynInd(T,vA,5);
        
        %RESIDUALS
        xB =x;
        if isfield(E,'Db');xB=bsxfun(@times,xB,dephaseBasis( E.Db.B, dynInd(E.Db.cr,vA(vA<=E.NMs),5), size(x)));end
        if isfield(E,'Dc');xB=bsxfun(@times,xB,dephaseRotation(dynInd(Tr,E.Dc.d,6),E.Dc.D));end   
        %xT=resPop(xB,6,[],4);        
        Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Tr,1,0,1,1);%Transform factors
        xT=sincRigidTransform(xB,Tf,1,E.Fof,E.Fob);        
        %xT=resPop(xT,4,[],6);
        
        E.vA=vA;
        E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end); %YB: set block sizes to the exact value as vA so that encode(x,E) computes 1 segment
        ry=encode(xT,E)-dynInd(y,indY,1);%Residuals
        EnPrev=EnPrev+gather(normm(ry,[],setdiff(1:5,dimR))); %YB: energy per coefficient and dimR (separable dimension)
        
        %JACOBIAN
        xT=dephaseRotationGradient(xB,dynInd(Tr,E.Dc.d,6),E.Dc.A); %YB: basis function must be before the encode! (or after decode as in MRes paper) 
        xT=resPop(xT,6:7,[],4);%YB: so that rotation does not apply to that dimension
        xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);
        xT=resPop(xT,4,ndS,6:7);% Place back in original dimensions
      
        %xT=dephaseRotationGradient(xT,dynInd(Tr,E.Dc.d,6),E.Dc.A); %YB:
        %reason why not here and before? see above
        
        J=encode(xT,E);
        dG=dG+gather(resPop(resPop(multDimSum(real(bsxfun(@times,J,conj(ry))),setdiff(1:5,dimR)),dimS,[],1),3,[],2));
        
        J=resPop(J,dimS,[],dimS(1)); %YB: is this reshaping in accordance with the ha creature? 
        for h=1:ndH %YB: hessian pairs in the Hessian matrix. Approximated as J'*J
            dH(h,:) = dH(h,:) + gather(resPop(multDimSum(real(bsxfun(@times,dynInd(J,ha(1,h),dimS(1)),conj(dynInd(J,ha(2,h),dimS(1))))),setdiff(1:5,dimR)),3,[],2));
        end
    end
    
    %UPDATE
    MH=single(eye(ndG));%YB: Matrix Hessian
    Dc.cupr=gather(Dc.c);Dc.cupr(:)=0;
    Dc.cup=Dc.c;Dc.cup(:)=0;
    Sff=E.Sf;
    while any(~flagw(:))   
        En=EnPrev;
        indNoConvFlag=find(~flagw(:));
        En(indNoConvFlag)=0; %YB: set energy to zero where no convergence     
       
        for a=indNoConvFlag' % %YB: going over the non-converged coefficients
            %BUILD HESSIAN MATRIX AND POTENTIAL UPDATE
            for h=1:ndH
                if ha(1,h)==ha(2,h);MH(ha(1,h),ha(2,h))=(1+E.Dc.w(a))*dH(h,a); % Only weight on diagonal
                else MH(ha(1,h),ha(2,h))=dH(h,a);MH(ha(2,h),ha(1,h))=dH(h,a);end
            end             
            GJ=single(double(MH)\double(dG(:,a)));
            upr=resPop(-(E.Dc.winit/E.Dc.w(a))*GJ,1,ndS,dimS);
            Dc.cupr=dynInd(Dc.cupr,a,3,upr);
        end
        if gpu;Dc.cupr=gpuArray(Dc.cupr);end
        Dc.cup=dynInd(Dc.cup,indNoConvFlag,3,dynInd(Dc.c,indNoConvFlag,3)+dynInd(Dc.cupr,indNoConvFlag,3));
        cup=dynInd(Dc.cup,indNoConvFlag,3);
        if NT(3)~=1
            yA=dynInd(y,indNoConvFlag,3);
            E.Sf=dynInd(Sff,indNoConvFlag,3);
        else
            yA=y;
        end
        Dup=sum(bsxfun(@times,cup,E.Dc.A),dimS(1)+1);%Residual phases per unit axial rotation
        
        %CHECK ENERGY REDUCTION
        for a=ET.oS(1):ET.bS(1):ET.dS(1);vA=a:min(a+ET.bS(1)-1,ET.dS(1));%vA=vA(vA<=E.NMs);
            indY=[];for c=1:length(vA);indY=horzcat(indY,E.nEc(vA(c))+1:E.nEc(vA(c)+1));end         
            Tr=dynInd(T,vA,5);
            
            %RESIDUALS
            xB =x;
            if isfield(E,'Db');xB=bsxfun(@times,xB,dephaseBasis( E.Db.B, dynInd(E.Db.cr,vA(vA<=E.NMs),5), size(x)));end
            if isfield(E,'Dc');xB=bsxfun(@times,xB,dephaseRotation(dynInd(Tr,E.Dc.d,6),Dup));end %YB: use the updated linear coefficients D
            %xT=resPop(xB,6,[],4);        
            Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Tr,1,0,1,1);%Transform factors
            xT=sincRigidTransform(xB,Tf,1,E.Fof,E.Fob);        
            %xT=resPop(xT,4,[],6);
            
            E.vA=vA;
            E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
            ry=encode(xT,E)-dynInd(yA,indY,1);%Residuals
            En=dynInd(En,indNoConvFlag,3,dynInd(En,indNoConvFlag,3)+gather(normm(ry,[],setdiff(1:5,dimR))));             
        end 
        En(E.Dc.w(:)>1e10)=EnPrev(E.Dc.w(:)>1e10);% YB: ? 
        flagw=En<=EnPrev;      

        if deb>=2;fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(EnPrev),sum(En));end
        E.Dc.w(indNoConvFlag(~flagw(indNoConvFlag)))=E.Dc.w(indNoConvFlag(~flagw(indNoConvFlag)))*multB;
        E.Dc.w(indNoConvFlag(flagw(indNoConvFlag)))=E.Dc.w(indNoConvFlag(flagw(indNoConvFlag)))/multA;                
        E.Dc.w(E.Dc.w<1e-8)=multA*E.Dc.w(E.Dc.w<1e-8);%To avoid numeric instabilities 

        if all(flagw(:))
            Dc.c=Dc.cup;
            Dc.R=sum(bsxfun(@times,Dc.c,E.Dc.A),dimS(1)+1)-Dc.D;
            Dc.D=Dc.D+Dc.R;
                        
            cMax=gather(multDimMax(abs(Dc.R),setdiff(1:dimS(end),dimR)));
            if NT(3)>1;E.Dc.cC(indRNoConv)=cMax<E.Dc.cL;else E.Dc.cC=cMax<E.Dc.cL;end%Convergence flag
            if deb>=2
                fprintf('Maximum change in phase per unit axial rotation (rad/rad): %s (tolerance: %s)\n',sprintf('%0.3f ',max(cMax)),sprintf('%0.3f ',E.Dc.cL));
                %fprintf('Number of converged groups of parameters: %d/%d\n',sum(1-E.Dc.cC),numel(E.Dc.cC));
            end        
        end
    end  
    if NT(3)>1        
        E.Dc.c=dynInd(E.Dc.c,indRNoConv,3,Dc.c);%Parameters of the DCT
        E.Dc.D=dynInd(E.Dc.D,indRNoConv,3,Dc.D);%Phases per unit axial rotation
    else
        E.Dc.c=Dc.c;
        E.Dc.D=Dc.D;
    end
    E.Tr=T;E.bS=bS;E.dS=dS;E.oS=oS;E.Sf=Sf;
    if gpu;E.Sf=gpuArray(E.Sf);end
    if isfield(E,'vA');E=rmfield(E,'vA');end
    Enmin=min(sum(EnPrev),sum(En));
    if ~isfield(E.Dc,'En');E.Dc.En=Enmin;else E.Dc.En=cat(1,E.Dc.En,Enmin);end%To keep a record of the achieved energy
end
