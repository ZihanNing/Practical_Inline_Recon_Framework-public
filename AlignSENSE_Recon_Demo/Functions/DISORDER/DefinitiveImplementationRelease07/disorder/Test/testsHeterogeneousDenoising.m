%addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06/Provisional/HeterogeneousSVDFilter'));
addpath(genpath('/home/lcg13/git/complexSVDShrinkageDWI/complexSVDShrinkageDWI'));

clearvars
%If gamma 0.85, r=1000,NC=2, s=0.00005*[1 2 5 10], k=100 OptShrink works
%better sometimes

%MATRIX DIMENSIONS
N=1024;
gamma=0.5;%0.5;%0.75;
M=floor(N*gamma);
MN=[M N];

%SIGNAL
s=10*sqrt(N);%(1:2:20);%0.5:0.5:10;%[100 100];%Column vector of size R, rows representing groups

%NC=2;
%r=0.001;
NC=1;
r=1;
s=bsxfun(@times,s,exp(linspace(log(r),0,NC)'));
s=flip(s,2);
[NC,R]=size(s);
dimHetero=1;
X=[];
isReal=2;

U=randU([M R],isReal);
V=randU([R N],isReal);
if dimHetero==2
    V=reshape(V,[R N/NC NC]);
    for n=1:NC;X=cat(2,X,matfun(@mtimes,bsxfun(@times,U,s(n,:)),V(:,:,n)));end
else
    U=reshape(U,[M/NC NC R]);
    for n=1:NC;X=cat(1,X,matfun(@mtimes,bsxfun(@times,squeeze(U(:,n,:)),s(n,:)),V));end
end
[W,WH]=buildWaveletMatrix(MN(dimHetero),'db4',4);
[F,FH]=buildStandardDFTM(MN(dimHetero),1);
W=W{1};WH=WH{1};
F=F{1};FH=FH{1};
if isReal==2;MS=eye(MN(dimHetero));
elseif isReal;MS=W';%We check a wavelet domain
else MS=F';%We check a Fourier domain
end
if dimHetero==1;X=MS'*X;else X=X*MS;end

%THIS WOULD BUILD AN OVERLAPPED WINDOW, IT HAS BEEN USELESS
%E=16;%2*N/NC;
%In=tukeywin(E,0);
%In=circconvmtx(In,N);
%subR=4;
%In=In(:,1:subR:end)';
%In=bsxfun(@times,In,1./sum(In,1));

%%%THERE ARE PROBLEMS COMPUTING THE AMSE---IT IS UNDERESTIMATING IT,
%%%PARTICULARLY WHEN ROTATING!!

%NOISE
k=[1 1];%[1 1000];
useMatrix=0;
C=linspace(k(1),k(2),M)';
if useMatrix;UG=randU([M M],isReal);L=UG'*diag(sqrt(C))*UG;else L=diag(sqrt(C));end
%C(1:M/2)=1;C(M/2+1:M)=k;%This is a way to make the weighting more appropriate!
G=plugNoise(X,isReal);

Ma=randi(2,[M 1],'like',real(X))-1;
%Ma(:)=1;
    
Xmax=max(max(X(Ma==1,:)))^2;
NMa=sum(Ma==1)*N;
%MATRIX AND PRECOMPUTATIONS
for ca=1:4
    if ca==1
        fprintf('Full signal\n');
    elseif ca==2
        fprintf('Masked signal\n');
        X=bsxfun(@times,Ma,X);     
    elseif ca==3
        fprintf('Masked signal and noise\n');
        X=bsxfun(@times,Ma,X);     
        Maf=Ma;
        Maf(Ma==0)=1e-6;%IT GENERALLY WORKS FOR 1E-2 SO THIS CAN BE FINE FOR SPATIAL MASKING        
        if useMatrix;L=diag(Maf)*UG'*diag(sqrt(C))*UG;else L=diag(Maf)*diag(sqrt(C));end
    elseif ca==4
        fprintf('Extracted signal\n');
        X=X(Ma==1,:);
        L=L(Ma==1,Ma==1);      
        G=G(Ma==1,:);
        MN=size(X);
        MS=MS(Ma==1,Ma==1);
        Ma=Ma(Ma==1);
    end
    Y=X+L*G;
    Y=double(Y);

    %TEST 1-WE DECORRELATE, ESTIMATE AND CORRELATE
    Ys=mldivide(L,Y);
    [amses,~,Xhats]=SVShrinkage(Ys);
    Xhats=L*Xhats;
    %amses=amses.*C;%May this give an estimate of the correlations in the denoised data?
    %amses=bsxfun(@times,amses,Ma);
    fprintf('MSE (dB): %.2f\n',10*log10(NMa*Xmax/multDimSum(bsxfun(@times,Ma,abs(Xhats-X)).^2,1:2)));
    %fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amses(:))));

    %TEST 2-WE APPLY OPTSHRINK
    %esd=ESDMixAndMix(C,N,[],2);
    esd=ESDSimulated(L'*L,N);
    [amse,~,Xhat]=SVShrinkage(Y,[],[],esd);
    %Xhat=bsxfun(@times,Ma,Xhat);
    %amse=bsxfun(@times,Ma,amse);
    fprintf('MSE (dB): %.2f\n',10*log10(NMa*Xmax/multDimSum(bsxfun(@times,Ma,abs(Xhat-X)).^2,1:2)));
    fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amse(:))));

    %TEST 3-WE APPLY STANDARD WEIGHTED ESTIMATION
    [Xhath,amseh]=heterogeneousSVShrinkage(Y,{{L(:)},[]});
    %Xhath=bsxfun(@times,Ma,Xhath);
    %amseh=bsxfun(@times,Ma,amseh);
    fprintf('MSE (dB): %.2f\n',10*log10(NMa*Xmax/multDimSum(bsxfun(@times,Ma,abs(Xhath-X)).^2,1:2)));
    fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amseh(:))));

    %TEST 4-WE APPLY LOCALIZED ESTIMATION
    NS=NC;
    InBase=ones(MN(dimHetero)/NS,1);
    In=[];
    for n=1:NS;In=blkdiag(In,InBase);end
    %In=cat(1,Ma'==1,Ma'==0);    
    
    if dimHetero==1;[Xhatll,amsell]=heterogeneousSVShrinkage(Y,{{L(:)},[]},{In,[]},{{MS(:)},[]});
    else [Xhatll,amsell]=heterogeneousSVShrinkage(Y,{{L(:)},[]},{[],In},{[],{MS(:)}});
    end
    %Xhatll=bsxfun(@times,Ma,Xhatll);
    %amsell=bsxfun(@times,Ma,amsell);
    fprintf('MSE (dB): %.2f\n',10*log10(NMa*Xmax/multDimSum(bsxfun(@times,Ma,abs(Xhatll-X)).^2,1:2)));
    fprintf('AMSE (dB): %.2f\n',10*log10(Xmax/mean(amsell(:))));
end