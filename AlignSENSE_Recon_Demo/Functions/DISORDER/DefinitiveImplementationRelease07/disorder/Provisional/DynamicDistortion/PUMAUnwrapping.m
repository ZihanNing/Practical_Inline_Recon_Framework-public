function B=PUMAUnwrapping(x,M,modal)

% PUMAUNWRAPPING calls the Phase Unwrapping MAx-flow/min-cut in [1] JM
% Bioucas-Dias, G Valadao. "Phase Unwrapping via Graph Cuts," IEEE Trans
% Image Proc, 16(3):698-709, March 2007 for dynamic field mapping
%   B=PUMAUNWRAPPING(X,M,MODAL)
%   * X is the complex data used to estimate the field
%   * W is a soft mask
%   * {MODAL} is the data modality (9->fMRI / 10->DWI), it defaults to 9
%   * B is the estimated field
%

if ~exist('modal','var') || isempty(modal);modal=9;end

%COMPILATION OF THE GC CODE---THIS SHOULD BE PRECOMPILED IN THE FUTURE
GCFold=fullfile(fileparts(mfilename('fullpath')),'../../PhaseUnwrapping','to_compile_mf2');
mex('-silent','-outdir',GCFold,fullfile(GCFold,'mf2.cpp'),fullfile(GCFold,'graph.cpp'),fullfile(GCFold,'maxflow.cpp'));

%UNWRAPPING PARAMETERS
N=size(x);N(end+1:5)=1;
ND=numDims(x);
lnorm=2;
w=[5;10;2;1;0.5];w=w(1:ND);    
%w=[5;10;1;0.1;0.5];w=w(1:ND);
%fprintf('Used weights:%s\n',sprintf(' %.2f',w));
c=eye(ND,ND); 
%TEMPORALLY AVERAGED UNWRAPPING
if N(4)>1
    xinit=mean(x,4);xinit=reshape(xinit,[N(1:3) N(5)]);
    winit=w;winit(4)=[];
    cinit=c(1:ND-1,1:ND-1);
    qinit=repmat(reshape(1-mean(M,4),[N(1:3) N(5)]),[ones(1,ND-1) ND-1]);
    zinit=angle(xinit);
    [~,k]=puma_hoND(zinit,lnorm,'c',cinit,'q',qinit,'w',winit);    
    k=reshape(k,[N(1:3) 1 N(5)]);
    k=repmat(k,[ones(1,3) N(4) 1]);
else
    k=zeros(N);
end

%4D UNWRAPPING
B=angle(x);
NB=size(B);
NQ=size(M);NQ(end+1:ND)=1;
radWin=1;    
if modal==10        
    q=repmat(1-M,[NB./NQ ND]);
    B=puma_hoND(B,lnorm,'c',c,'q',q,'w',w,'k0',k);
else
    Baux=B;                        
    if radWin==0
        w(4)=[];c=c(1:ND-1,1:ND-1);
        q=repmat(reshape(1-mean(M,4),[N(1:3) N(5)]),[ones(1,ND-1) ND-1]);
    else                                                
        q=repmat(1-M,[NB./NQ ND]);
    end        
    for n=1:N(4)
        vW=n-radWin:n+radWin;vW=mod(vW-1,N(4))+1;
        B=dynInd(B,n,4,dynInd(puma_hoND(dynInd(Baux,vW,4),lnorm,'c',c,'q',dynInd(q,vW,4),'w',w),radWin+1,4));
    end
end