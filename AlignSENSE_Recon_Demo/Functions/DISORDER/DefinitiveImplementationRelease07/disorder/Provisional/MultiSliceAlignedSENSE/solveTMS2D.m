function [T,w,flagw,outlD]=solveTMS2D(x,y,T,H,S,Ak,xkGrid,kkGrid,kGrid,nT,w,flagw,NzSlab,outlP,thplc,BlSz,winit)

%SOLVETMS2D estimates motion during acquisition
%   [T,W,FLAGW,OUTLD]=SOLVETMS2D(X,Y,T,H,S,AK,XKGRID,KKGRID,KGRID,NT,W,FLAGW,NZSLAB,GPU,OUTLP,THPLC,BLSZ)
%   computes the best T for a given x.
%   * X is the reconstructed image
%   * Y is the measured data
%   * T are the input transform parameters
%   * H is the slice profile filter
%   * S are the sensitivities
%   * AK is a sampling mask in the phase encoding direction
%   * XKGRID is a grid of points in the spatial-spectral domain
%   * KKGRID is a grid of points in the spectral-spectral domain
%   * KGRID is a grid of points in the spectral domain
%   * NT is the number of iterations of the Hessian method (please note 
%   that this code has not been tested for NT>1)
%   * W are the weights of the Newton's method
%   * FLAGW is a flag to indicate whether to increment or decrement the
%   weights
%   * NZSLAB is the number of slices of the slab
%   * OUTLD is the threshold for shot rejection
%   * THPLC is a flag that determines whether to estimate through-plane
%   motion terms
%   * BLSZ is the size of blocks for gpu processing
%   ** T is the output transform parameters
%   ** W are the weights of the Newton's method
%   ** FLAGW is a flag to indicate whether to increment or decrement the
%   weights
%   ** OUTLD is a mask for shot rejection
%

gpu=isa(x,'gpuArray');
NS=size(S);NY=size(y);

multA=1.2;%Factor to divide the weight (that regularizes the Hessian matrix) when H(end)<H(end-1)
multB=2;%Factor to multiplicate the weight (that regularizes the Hessian matrix) when H(end)>=H(end-1)

NT=size(T);
Tup=T;

a(1,:)=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6];a(2,:)=[1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];

NHe=size(a,2);

Eprev=single(zeros(NT(1:6)));E=single(zeros(NT(1:6)));
dHe=single(zeros([NHe NT(5:6)]));
dH=single(zeros([NT(7) NT(5:6)]));

if gpu;Eprev=gpuArray(Eprev);E=gpuArray(E);end

%Masking (usually aids convergence although it may be disabled)
MaskThreshold=0.2;
xabs=abs(x);
meanBody=mean(xabs(:));
xabs=xabs>meanBody*MaskThreshold;
x=x.*xabs;

NZR=size(x);
NZR(end+1:5)=1;
yZR{2}{2}=zeros([NZR(1:2) NzSlab NZR(4:5) NZR(3)],'like',x);
yZR{1}{2}=zeros([NZR(1:4) NT(5)],'like',x); 
x=extractSlabsOld(x,NzSlab,1,1,yZR); 

%Iterations
for n=1:nT
    %Update the weights                   
    w(flagw==2)=w(flagw==2)/multA;
    w(flagw==1)=w(flagw==1)*multB;
    for s=1:BlSz:NT(5);vS=s:min(s:s+BlSz-1,NT(5));
        [et,etg,eth]=precomputeFactors3DTransform(xkGrid,kkGrid,kGrid,dynInd(T,vS,5),1,2); 
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                etg{1}{m}=gpuArray(etg{1}{m});
                for l=2:3
                    et{l}{m}=gpuArray(et{l}{m});etg{l}{m}=gpuArray(etg{l}{m});eth{l}{m}=gpuArray(eth{l}{m});
                end
            end
            for m=1:6;eth{1}{m}=gpuArray(eth{1}{m});end
        end        
        if s==1
            NX=size(x);
            [F,FH]=buildStandardDFTM(NX,0,gpu);
        end        
        [xT,xB]=transform3DSinc(x,et,1,F,FH);
        xT=encoding(xT);    
        xT=bsxfun(@minus,xT,y);
        xT=bsxfun(@times,xT,dynInd(Ak,vS,5));
        xTH=conj(xT);    
        Eprev=dynInd(Eprev,vS,5,permute(multDimSum(real(xT.*xTH),[1 2 4 6]),[1 2 6 4 5 3]));     

        [G,GB,GC]=transform3DSincGradient(xB,et,etg,0,F,FH);
        for m=1:NT(7)                   
            G{m}=encoding(G{m}); 
            G{m}=bsxfun(@times,G{m},dynInd(Ak,vS,5));
            GH{m}=conj(G{m});
        end
        for m=1:NHe%We don't use the Hessian!         
            GG=transform3DSincHessian(xB,GB,GC,et,etg,eth,m,F,FH);
            GG=encoding(GG); 
            GG=bsxfun(@times,GG,dynInd(Ak,vS,5));
            GG=real(GG.*xTH);
            GG=GG+real(G{a(1,m)}.*GH{a(2,m)});  
            %GG=real(G{a(1,m)}.*GH{a(2,m)});
            dHe(m,vS,:)=gather(permute(multDimSum(GG,[1 2 4 6]),[1 5 3 2 4]));
        end
        
        for m=1:NT(7)             
            G{m}=real(bsxfun(@times,G{m},xTH)); 
            dH(m,vS,:)=gather(permute(multDimSum(G{m},[1 2 4 6]),[1 5 3 2 4]));  
        end
    end
   MHe=single(1000*eye(NT(7)));      
    for t=1:NT(6)
        for s=1:NT(5)
            for k=1:NHe
                if a(1,k)==a(2,k)
                    MHe(a(1,k),a(2,k))=dHe(k,s,t)+w(1,1,1,1,s,t);                            
                else
                    MHe(a(1,k),a(2,k))=dHe(k,s,t);
                    MHe(a(2,k),a(1,k))=dHe(k,s,t);
                end              
            end                        
            dH(:,s,t)=single(double(MHe)\double(dH(:,s,t)));
        end
    end

%     MHe=single(eye(NT(7)));
%     for t=1:NT(6)
%         for s=1:NT(5)
%             for k=1:NHe
%                 if a(1,k)==a(2,k)
%                     MHe(a(1,k),a(2,k))=dHe(k,s,t)*(1+w(1,1,1,1,s,t));                            
%                 else
%                     MHe(a(1,k),a(2,k))=dHe(k,s,t);
%                     MHe(a(2,k),a(1,k))=dHe(k,s,t);
%                 end              
%             end                        
%             dH(:,s,t)=(winit/w(1,1,1,1,s,t))*single(double(MHe)\double(dH(:,s,t)));
%         end
%     end
    
    
    if thplc==0;dH([1 2 4],:,:)=0;end
    if thplc==1;dH([3 5 6],:,:)=0;end    
    
    Tup=permute(dH,[4 5 6 7 2 3 1]);     
    Tup=T-Tup;            
    
    for s=1:BlSz:NT(5);vS=s:min(s:s+BlSz-1,NT(5));
        et=precomputeFactors3DTransform(xkGrid,kkGrid,kGrid,dynInd(Tup,vS,5),1,0);
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                for l=2:3;et{l}{m}=gpuArray(et{l}{m});end
            end
        end
        xT=transform3DSinc(x,et,1,F,FH);          
        xT=encoding(xT);        
        xT=bsxfun(@minus,xT,y);
        xT=bsxfun(@times,xT,dynInd(Ak,vS,5));
        xTH=conj(xT);
        E=dynInd(E,vS,5,permute(multDimSum(real(xT.*xTH),[1 2 4 6]),[1 2 6 4 5 3]));  
    end
    flagw(E<Eprev)=2;    
    flagw(E>=Eprev)=1;
    for s=1:NT(7)
        TauxA=dynInd(T,s,7);
        TauxB=dynInd(Tup,s,7);
        TauxA(flagw==2)=TauxB(flagw==2);
        T=dynInd(T,s,7,TauxA);
    end
    Eaux=min(E,Eprev);
    Eaux=2*Eaux./(Eaux(:,:,:,:,:,[2:end 1])+Eaux(:,:,:,:,:,[end 1:end-1]));
    outlD=permute((Eaux<outlP),[1 2 6 3 5 4]);
end

function x=encoding(x)
    x=filtering(x,H); 
    x=extractSlabsOld(x,NzSlab,0,1,yZR);  
    x=bsxfun(@times,x,S);
    %x=fold(x,1,NS(1),NY(1));
    x=fftGPU(x,1,[],NY(1));
end

end
