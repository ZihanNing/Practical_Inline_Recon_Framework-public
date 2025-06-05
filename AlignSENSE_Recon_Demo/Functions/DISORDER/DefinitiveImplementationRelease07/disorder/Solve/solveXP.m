function [x,P]=solveXP(x,y,E,EH,A,C,R,nIt,tolType,tol,deb)

%SOLVEXP   Performs a simultaneous estimation of the phase and the reconstructions
%   [X,P]=SOLVEXP(X,Y,E,EH,{A},{C},{R},{NIT},{TOLTYPE},{TOL},{DEB})
%   X is the initial data reconstruction
%   Y is the array of measures
%   E is an encoding structure
%   EH is a decoding structure
%   {A} is a preconditioner structure
%   {C} is a constrain structure
%   {R} is a regularizer structure
%   {NIT} is the maximum number of iterations
%   {TOLTYPE} is the type of tolerance used for convergence
%   {TOL} is the tolerance set for convergence
%   {DEB} indicates whether to print information about convergence
%   X is the reconstruction
%   P are the estimated phase parameters
%
%The NormwiseBackwardError stopping condition is based on Petr Tichy. 
%ON ERROR ESTIMATION IN THE CONJUGATE GRADIENT METHOD: NORMWISE BACKWARD 
%ERROR, Proceedings of ALGORITMY, pp. 323-332 2016.

if nargin<5;A=[];end
if nargin<6;C=[];end
if nargin<7;R=[];end
if nargin<8 || isempty(nIt);nIt=300;end%300;end
if nargin<9 || isempty(tolType);tolType='Energy';end%tolType='NormwiseBackward2Error'
if nargin<10 || isempty(tol);tol=1e-5;end%5e-3;end%1e-2;end
if nargin<11 || isempty(deb);deb=1;end

gpu=isa(x,'gpuArray');if gpu;gpuF=2;else gpuF=0;end

%HARDCODED PARAMETERS
nExtern=100;%Maximum number of external iterations
EP.tolP=2e-3;%0.01;%Tolerance for P
EP.tolEn=1e-4;%Tolerance in terms of the energy-We see no differences below this threshold. Guess a problem would be when the decrease is too low by chance...
multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)

%WE START BY EXTRACTING THE PARAMETERS OF THE PHASE MAP
[~,P]=ridgeDetection(E.Gh.Pf,setdiff(1:2,E.pe),16);
yF=y;EF=E;EHF=EH;AF=A;CF=C;RF=R;PF=P;xF=x;

MB=E.Uf{3}.NX/E.Uf{3}.NY;%Multiband factor
NP=size(P);%Constant and linear phase/Odd and even lines/Number of slices
EP.dimP=[1 NP(E.pe) E.Uf{3}.NY];
EP.cXsel=false([1 1 E.Uf{3}.NY]);
EP.cP=false(EP.dimP);%Convergence of estimates
EP.cX=false(EP.dimP);%Convergence of reconstruction
EP.flagw=single(zeros(EP.dimP));%Convergence of reconstruction
ndG=NP(setdiff(1:2,E.pe))*MB;
ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
ndH=size(ha,2);
EP.EnPrev=single(zeros(EP.dimP));EP.En=EP.EnPrev;
EP.winic=single(0.01);%Weights of LM
EP.w=EP.winic*single(ones(EP.dimP));%Weights of LM
EPF=EP;
for n=1:nExtern    
    %INITIALIZE                
    EPF.cXsel=all(EPF.cX,2);
    sl=find(~EPF.cXsel(:))';

    [y,P,x,E,EH,A,C,R,EP]=extractExcitations(yF,PF,xF,EF,EHF,AF,CF,RF,EPF,sl);     
    EP.flagw(:)=0;EP.cP(:)=0;
    
    NP=size(P);%Constant and linear phase/Odd and even lines/Number of slices
    dP=single(zeros([NP(setdiff(1:2,E.pe))*MB NP(E.pe) E.Uf{3}.NY]));
    ddP=single(zeros([ndH NP(E.pe) E.Uf{3}.NY]));
    NS=size(E.Sf);NS(end+1:4)=1;%Size of the coils   
    if E.pe==1;NJ=[E.Uf{1}.NY NS(2) E.Uf{3}.NY NS(4) NP(setdiff(1:2,E.pe))*MB];else NJ=[NS(1) E.Uf{2}.NY E.Uf{3}.NY NS(4) NP(setdiff(1:2,E.pe))*MB];end  
    J=zeros(NJ,'like',x);%Jacobian
    Pprev=P;
    
    %COMPUTE ENERGY
    E.Gh.Pf=phaseMap(P,setdiff(1:2,E.pe),size(E.Gh.Pf));   
    [E,EH]=computeEGhA(E,EH);

    xE=conj(encode(x,E)-y);
    dS=multDimSum(real(xE.*conj(xE)),[setdiff(1:2,E.pe) 4]);
    for d=1:length(E.Gh.nD);e=E.Gh.nD(d);
        indPE=(E.Gh.Mf==e);        
        EP.EnPrev=dynInd(EP.EnPrev,d,2,gather(sum(dynInd(dS,indPE,E.pe),E.pe)));
    end
    
    %GRADIENT    
    Pfg=phaseMapGradient(P,setdiff(1:2,E.pe),size(E.Gh.Pf));
    %JACOBIAN
    c=1;    
    for m=1:MB
        xJ=x;NX=size(x);NX(end+1:6)=1;
        xJ=reshape(xJ,[NS(1:2) E.Uf{3}.NY MB 1 NX(6)]);
        xJ=dynInd(xJ,setdiff(1:MB,m),4,0);
        xJ=reshape(xJ,[NS(1:3) 1 1 NX(6)]);
        for l=1:NP(setdiff(1:2,E.pe))
            E.Gh.Pf=Pfg{l};
            [E,EH]=computeEGhA(E,EH);
            J=dynInd(J,c,5,encode(xJ,E));
            c=c+1;
        end
    end

    %DERIVATIVES
    for g=1:ndG
        dS=multDimSum(real(dynInd(J,g,5).*xE),[setdiff(1:2,E.pe) 4]);
        for d=1:length(E.Gh.nD);e=E.Gh.nD(d);
            indPE=(E.Gh.Mf==e);
            dP=dynInd(dP,[g d],1:2,gather(sum(dynInd(dS,indPE,E.pe),E.pe)));
        end
    end    
    %DERIVATIVE OUTER PRODUCT    
    for h=1:ndH
        dS=multDimSum(real(dynInd(J,ha(1,h),5).*conj(dynInd(J,ha(2,h),5))),[setdiff(1:2,E.pe) 4]);
        for d=1:length(E.Gh.nD);e=E.Gh.nD(d);
            indPE=(E.Gh.Mf==e);
            ddP=dynInd(ddP,[h d],1:2,gather(sum(dynInd(dS,indPE,E.pe),E.pe)));
        end
    end
    dPA=reshape(dP,[ndG NP(E.pe)*E.Uf{3}.NY]);
    ddPA=reshape(ddP,[ndH NP(E.pe)*E.Uf{3}.NY]);
                  
    fina=0;
    MH=single(eye(ndG));    
    while ~fina
        dPEff=dP;
        dPEff(:)=0;
        for a=find(~EP.cP(:) & EP.w(:)<1e10)'           
            for h=1:ndH
                if ha(1,h)==ha(2,h)
                    MH(ha(1,h),ha(2,h))=(1+EP.w(a))*ddPA(h,a);
                else
                    MH(ha(1,h),ha(2,h))=ddPA(h,a);
                    MH(ha(2,h),ha(1,h))=ddPA(h,a);
                end
            end 
            MH=MH+mean(diag(MH))*single(eye(ndG))/1e6;%To avoid numerical instabilities
            dPEff(:,a)=-(EP.winic/EP.w(a))*single(double(MH)\double(dPA(:,a)));
        end       
        if gpu;dPEff=gpuArray(dPEff);end
        dPEff=reshape(dPEff,[NP(setdiff(1:2,E.pe)) MB NP(E.pe) E.Uf{3}.NY]);
        if E.pe==1;dPEff=permute(dPEff,[3 1 4 2]);else dPEff=permute(dPEff,[1 3 4 2]);end        
        dPEff=reshape(dPEff,NP);        
        Pup=P+dPEff;
                
        %CHECK ENERGY REDUCTION
        E.Gh.Pf=phaseMap(Pup,setdiff(1:2,E.pe),size(E.Gh.Pf));        
        [E,EH]=computeEGhA(E,EH);
        xE=conj(encode(x,E)-y);  
        dS=multDimSum(real(xE.*conj(xE)),[setdiff(1:2,E.pe) 4]);        
        for d=1:length(E.Gh.nD);e=E.Gh.nD(d);
            indPE=(E.Gh.Mf==e);
            EP.En=dynInd(EP.En,d,2,gather(sum(dynInd(dS,indPE,E.pe),E.pe)));
        end
        EP.En(EP.w(:)>1e4)=EP.EnPrev(EP.w(:)>1e4);%THIS HAS BEEN CHANGED TO ACCELERATE
        %EP.En(EP.w(:)>1e10)=EP.EnPrev(EP.w(:)>1e10);
        EP.En(EP.cP(:))=EP.EnPrev(EP.cP(:));
        EP.flagw(EP.En<=EP.EnPrev)=1;EP.flagw(EP.cP)=2;
        if deb>=2;fprintf('Iter %d - Energy before: %0.6g / Energy after: %0.6g\n',n,sum(EP.EnPrev(:)),sum(EP.En(:)));end

        if any(~EP.flagw(:)) 
            EP.w(EP.flagw==0)=EP.w(EP.flagw==0)*multB;
        else
            EP.w(EP.flagw==1)=EP.w(EP.flagw==1)/multA;
            EP.w(EP.w<1e-8)=multA*EP.w(EP.w<1e-8);%To avoid numeric instabilities 
            EP.cP(EP.flagw==1)=1;
            [P,Pup]=parUnaFun({P,Pup},@reshape,[NP(1:2) E.Uf{3}.NY MB]);
            [P,Pup]=parUnaFun({P,Pup},@permute,[setdiff(1:2,E.pe) 4 E.pe 3]);
            [P,Pup]=parUnaFun({P,Pup},@reshape,[NP(setdiff(1:2,E.pe)) MB NP(E.pe)*E.Uf{3}.NY]);
            
            P=dynInd(P,EP.flagw==1,3,dynInd(Pup,EP.flagw==1,3));
            P=reshape(P,[NP(setdiff(1:2,E.pe)) MB NP(E.pe) E.Uf{3}.NY]);
            if E.pe==1;P=permute(P,[3 1 4 2]);else P=permute(P,[1 3 4 2]);end
            P=reshape(P,NP);       
        end
        if all(EP.flagw(:));fina=1;end
    end    

    %CHECK CONVERGENCE
    Pup=P-Pprev;
    Pup=dynInd(Pup,2,setdiff(1:2,E.pe),dynInd(Pup,2,setdiff(1:2,E.pe))*NS(setdiff(1:2,E.pe)));
    Pup=reshape(Pup,[NP(setdiff(1:2,E.pe))*MB NP(E.pe) E.Uf{3}.NY]);
    Pup=permute(Pup,[2 3 1]);
    Pup=max(abs(Pup),[],3);
    EP.cX=gather(Pup)<pi*EP.tolP | shiftdim(gather((EP.EnPrev-EP.En)./EP.EnPrev),1)<EP.tolEn;
    
    %SOLVE FOR X
    E.Gh.Pf=phaseMap(P,setdiff(1:2,E.pe),size(E.Gh.Pf));EH.Gh.Pb=conj(E.Gh.Pf);
    [E,EH]=computeEGhA(E,EH);
    x=CGsolver(y,E,EH,A,C,R,x,nIt,tolType,tol,[],[],deb);
    %ASSIGN ESTIMATES
    [PF,xF,EF,EHF,EPF]=assignExcitations(PF,xF,EF,EHF,EPF,P,x,E,EH,EP,sl);    

    if all(EPF.cX(:));break;end        
end
if n>=nExtern;fprintf('Maximum number of iterations reached in Nyquist ghosting estimation\n');end
x=xF;
P=PF;