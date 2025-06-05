function [B,W]=CNCGUnwrappingOld(x,lp,voxSiz,upsampling,weightFunc)

%CNCGUNWRAPPING calls an lp-norm phase unwrapping (Convex Norm Conjugate
%Gradient) based on [1] DC Ghiglia, LA Romero. "Minimum Lp-norm 
%two-dimensional phase  unwrapping," J Opt Soc Am A, 13(10):1999-2013, Oct 
%1999, which is pre-weighted using the magnitude data as proposed in [2] DC 
%Ghiglia and LA Romero. "Robust two-dimensional weighted and unweighted
%phase unwrapping that uses fast transforms and iterative methods," J Opt
%Soc Am A, 11(1):107-117, Jan 1994, for dynamic field mapping
%   [B,W]=CNCGUNWRAPPING(X,{LP},{VOXSIZ},{UPSAMPLING})
%   * X is the complex data used to estimate the field
%   * {LP} is the norm of the problem
%   * {VOXSIZ} are the voxel sizes of the data
%   * {UPSAMPLING} upsampling factor for unwrapping. Defaults to 1
%   * B is the estimation of the unwrapped field
%   * W is the distance between the unwrapped and the wrapped field
%

NXOr=size(x);
x=squeeze(x);%%%CHANGE TO SAVE MEMORY
NXOrb=size(x);
if nargin>=3;voxSiz(NXOrb==1)=[];end
if nargin<2 || isempty(lp);lp=1;end
ND=numDims(x);
if nargin<3 || isempty(voxSiz);voxSiz=ones(1,ND);end
if nargin<4 || isempty(upsampling);upsampling=1;end
if nargin<5;weightFunc='Magnitude';end

gpu=isa(x,'gpuArray');

%UNWRAPPING PARAMETERS
N=size(x);
w=voxSiz;w(end+1:ND)=1;
%fprintf('Used weights:%s\n',sprintf(' %.2f',w));

%upSample=[2.5 2.5 1 1];
%upSample=ones(1,ND);
upSample=[upsampling upsampling ones(1,ND-2)];upSample=upSample(1:ND);
w=w./upSample;

%PHASE UNWRAPPING ENCODING STRUCTURES
%Encoder/Decoder

Nup=round(N.*upSample);
x=resampling(x,Nup,[],2*ones(1,ND));
%Precondition
A.mi=ones(1,ND);
A.Ps=-buildFilter(Nup*2,'2ndFiniteDiscrete',w,0,0,A.mi);
if gpu;A.Ps=gpuArray(A.Ps);end
A.Ps=real(A.Ps);
A.Ps(1)=1e9;%For inversion
A.Ps=A.Ps.^(-1);

%PHASE ENCODING PRECOMPUTATIONS
B=angle(x);Babs=abs(x);x=[];

for n=1:ND
    NN=size(Babs,n);
    padEl=zeros(1,ND+1);padEl(n)=1;
    Baux=diff(B,1,n);
    Baux=padarray(wrapToPi(Baux)/w(n),padEl,0,'post');
    Bauxabs=ones(1,'like',Babs);
    if ~isempty(strfind(weightFunc,'Magnitude'));Bauxabs=bsxfun(@times,Bauxabs,sqrt(Babs.*dynInd(Babs,[2:NN 1],n)));end     
    if ~isempty(strfind(weightFunc,'Gradient'));Bauxabs=bsxfun(@times,Bauxabs,abs(cos(angle(exp(1i*B).*exp(-1i*dynInd(B,[2:NN 1],n))))));end      
    %if ~isempty(strfind(weightFunc,'Gradient'));Bauxabs=bsxfun(@times,Bauxabs,abs(cos(angle(Babs.*exp(1i*B)-dynInd(Babs.*exp(1i*B),[2:NN 1],n)))));end     
    
    if n==1
        Bn=Baux;
        EH.Mb=Bauxabs;
    else
        Bn=cat(ND+1,Bn,Baux);
        EH.Mb=cat(ND+1,EH.Mb,Bauxabs);
    end   
end
if isfield(EH,'Mb');lp=2;end

E.Ff.Fd=1;EH.Fb.Fd=1;
E.Ff.w=w;EH.Fb.w=w;
E.Ff.ND=ND;EH.Fb.ND=ND;
E.Sm=1;

Bor=B;
B=Bn;Bn=[];Baux=[];Bauxabs=[];
if E.Sm;B=gather(B);end
%C.mV=1:ND;
C=[];
R=[];


%PHASE UNWRAPPING SOLVER
Bor=gather(Bor);
if lp==2;B=CGsolver(B,E,EH,A,C,R);%%%THIS HAS BEEN RECENTLY ADDED BECAUSE EH.MB DOES NOT WORK WHEN USING IRWLSSOLVER
else B=IRWLSsolver(B,E,EH,lp,1:ND+1,[],A,C,R);
end;A=[];
if gpu;Bor=gpuArray(Bor);end
W=B-Bor;Bor=[];
[B,W]=parUnaFun({B,W},@resampling,N);
[B,W]=parUnaFun({B,W},@real);
[B,W]=parUnaFun({B,W},@reshape,NXOr);



