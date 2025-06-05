function [xou, xou_notsum] =decode(x,EH,E, notsum)

%DECODE   Applies a given decoding to the data
%   XOU=DECODE(X,EH,{E})
%   * X is the data before decoding
%   * EH is the decoding structure
%   * {E} is the encoding structure (memory demanding or replicated 
%   information may be stored in this structure)
%   ** XOU is the data after decoding
%

if nargin < 4 || isempty (notsum) ; notsum = 0;end %YB for in GDsolver
xou_notsum = [];
if notsum ==2; E = rmfield(E,'Tr');end%avoid computing inverse transformations twice

if isempty(EH);xou=x;return;end
if nargin<3;E=[];end
inBlSz=1e3;

if ~isfield(E,'bS')%Block sizes by default
    E.bS=[size(x,5) size(x,4)];E.dS=[size(x,5) size(x,4)];E.oS=[1 1];
    if isfield(E,'Sf');E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end
    if isfield(E,'Tr');E.bS(1)=E.NSe;E.dS(1)=E.NSe;end
end
if isfield(E,'Zf') && isempty(E.Zf) && isfield(E,'Sf');E.oS(2)=1;E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end

if isfield(EH,'Mb');x=bsxfun(@times,x,EH.Mb);end%Mask / Weighting
xou=x;
for a=E.oS(1):E.bS(1):E.dS(1);vA=a:min(a+E.bS(1)-1,E.dS(1));
    xT=x;
    xouT=xT;
    for b=E.oS(2):E.bS(2):E.dS(2);vB=b:min(b+E.bS(2)-1,E.dS(2));
        xS=dynInd(xT,vB,4);

        xouS=xS;
        
        %SLICE MASK
        if isfield(E,'Bm') && ~isempty(E.Bm);xouS=bsxfun(@times,xouS,E.Bm);end
        
        %FOURIER DOMAIN (MULTISLICE)
        if isfield(EH,'Fms') && ~isempty(EH.Fms)
            if isfield(EH,'We');xouS=bsxfun(@times,dynInd(EH.We,vA,5),dynInd(xouS,vA,5));else xouS=dynInd(xouS,vA,5);end
            xouS=matfun(@mtimes,dynInd(EH.Fms,vA,5),xouS);
        end
        %YB: At this point, x is the k-space and the first 2 demensions are reshaped to a 1D signal
        
        %FOURIER DOMAIN (SEGMENTS) 
        if isfield(E,'Fs') && ~isempty(E.Fs)
            for c=1:length(vA)
                xR=dynInd(xS,E.nEc(vA(c))+1:E.nEc(vA(c)+1),1);%YB: Now the indiced E.nEc are important -  starting indexes of the Fourier encoding matrices
                if isfield(EH,'We') && ~isempty(EH.We);xR=bsxfun(@times,xR,EH.We(E.nEc(vA(c))+1:E.nEc(vA(c)+1)));end %YB: We only in the decode step!
                xouR=[];
                if ~isempty(E.Fs{1}{vA(c)})
                    Nec=size(E.Fs{1}{vA(c)},1);
                    for d=1:inBlSz:Nec;vD=d:min(d+inBlSz-1,Nec);   %YB: x still 1D list of k-space e.g. 15x1                 
                        if size(E.Fs,1)==2;xQ=bsxfun(@times,dynInd(xR,vD,1),conj(dynInd(E.Fs{2}{vA(c)},vD,1)));end  %YB: by multiplying e.g. 15x70 and 
                        xQ=aplGPU(dynInd(E.Fs{1}{vA(c)},vD,1)',xQ,1); %YB: Matrix multiplication applies ifft in first direction and thus 60x7 = image dimension
                        if d==1;xouR=xQ;else xouR=xouR+xQ;end
                    end
                end;xQ=[];
                if c==1 || isempty(xouS);xouS=xouR;elseif ~isempty(xouR);xouS=cat(5,xouS,xouR);end %YB: Stack different motion states in 5th dim over a iterations
            end;xR=[];xouR=[];
        end
        
        if isfield(E,'Es') && E.Es==1 && isfield(E,'Ef');xouS=aplGPU(E.Ef,xouS,E.pe);end%Fourier
        
        %WAVE ENCODING - would come here
        if isfield(E,'Wave')
            xouS = fftshiftGPU ( fft(xouS,[],E.Wave.RO), E.Wave.RO); %To hybrid space - fft in readout
            xouS = bsxfun(@times, xouS, conj(E.Wave.PSF));%Apply Wave psf
            xouS = ifft ( ifftshiftGPU(xouS, E.Wave.RO),[] ,E.Wave.RO); %Back to image space - ifft in readout
            xouS = resampling( xouS, E.Wave.N, 2);%remove over-sampling
        end
        
        %SENSE
        if isfield(EH,'Ub')%Sense unfolding (third dimension)             
            for n=3
                if isfield(EH,'Bb');xouS=repmat(xouS,[1 1 EH.Ub{n}.NX/EH.Ub{n}.NY]);elseif ~isempty(EH.Ub{n});xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
            end   
        end
        
        %MB/NYQUIST
        if isfield(EH,'Gh');xouS=bsxfun(@times,xouS,EH.Gh.Ab);end%Nyquist ghosting decoding
        if isfield(EH,'Bb');xouS=bsxfun(@times,xouS,dynInd(EH.Bb,E.cc,6));end%MB decoding
        if isfield(EH,'Eb');xouS=aplGPU(EH.Eb,xouS,E.pe);end%Space
        
        if isfield(EH,'Ub')%Sense unfolding (first two dimensions)
            for n=2:-1:1                
                if ~isempty(EH.Ub{n})
                    if EH.Ub{n}.NX>EH.Ub{n}.NY;xouS=ifold(xouS,n,EH.Ub{n}.NX,EH.Ub{n}.NY,EH.UAb{n});end
                    if EH.Ub{n}.NX<EH.Ub{n}.NY
                        NXOUS=size(xouS);NXOUS(n)=EH.Ub{n}.NX;
                        xouS=resampling(xouS,NXOUS,2);
                    end
                end
            end
        end
        %COIL PROFILES
        if isfield(E,'Sf')
            if isfield(E,'coilLoading') && isfield(E.coilLoading,'coilIdx') && ~isempty(E.coilLoading.coilIdx)
                Saux = dynInd(E.Sf,{vB, E.coilLoading.coilIdx(vA(vA<=E.NMs))},4,5);
            elseif isfield(E,'coilLoading') && isfield(E.coilLoading,'f') && isfield(E,'Tr') && any(E.Tr(:)~=0)
                Saux = bsxfun(@times,dynInd(E.Sf,vB,4), coilLoading(dynInd(E.Tr,vA(vA<=E.NMs),5), vB, E.coilLoading.f));
            else
                Saux=dynInd(E.Sf,vB,4);
            end
            if isa(xouS,'gpuArray');Saux=gpuArray(Saux);end
            xouS=bsxfun(@times,xouS,conj(Saux));
            if isfield(E,'Zf') && ~isempty(E.Zf);xouS=bsxfun(@times,xouS,dynInd(E.Zf,vB,4));end
            if ~isfield(E,'Zf') || ~isempty(E.Zf);xouS=sum(xouS,4);end
        end%Sensitivities
        if b==E.oS(2) || isempty(xouT);xouT=xouS;elseif isfield(E,'Sf') && ~isempty(xouS);xouT=xouT+xouS;elseif ~isempty(xouS);xouT=cat(4,xouT,xouS);end %YB: need to sum over channels!
    end;xS=[];xouS=[];
    
    %SLAB EXTRACTION
    if isfield(E,'ZSl');xouT=extractSlabs(xouT,abs(E.ZSl),1,0);end
    
    %FILTERING (GENERALLY FOR SLICE RECOVERY)
    if isfield(E,'Sp');xouT=filtering(xouT,conj(E.Sp));end
    
    %B0 SHIM MODULATION
    if isfield(E,'Shim');xouT=bsxfun(@times,xouT,conj(exp(+1i*2*pi*E.Shim.TE *E.Shim.B0)));end   
        
    %RIGID TRANSFORM (MOTION STATES)
    if  isfield(E,'NMs') && any(vA<=E.NMs)   %YB reversed if-conditions
        xinT=dynInd(xouT,vA<=E.NMs,5);

        %DEPHASING (in scanner reference frame)
        if isfield(E,'Dbs');xinT=bsxfun(@times,xinT,conj(dephaseBasis(E.Dbs.B, dynInd(E.Dbs.cr,vA(vA<=E.NMs),5), size(xinT),E.Dbs.TE)));end
            
        %B1 modulation
        if isfield(E,'B1m');xinT=bsxfun(@times,xinT,exp(dephaseBasis(E.B1m.B, dynInd(E.B1m.cr,vA(vA<=E.NMs),5), size(xinT),[],1)));end   
    
        if notsum ==2 &&  ~isempty(xinT) % YB: Used in GDsolverScannerRef.m and GDsolverB1.m
                if isempty(xou_notsum); xou_notsum = xinT;
                else xou_notsum=cat(5,xou_notsum,xinT);end 
        end
        if isfield(E,'Tr') && ~isempty(E.Tr)   

            if any(E.Tr(:)~=0)           
                Tr=dynInd(E.Tr,vA(vA<=E.NMs),5);                                  
                %TRANSFORM
                if isfield(EH,'Tb');Tb=extractFactorsSincRigidTransform(EH.Tb,vA(vA<=E.NMs),5);
                else; Tb=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Tr,0,0,1,1); %YB: Take the inverse transform
                end
                
                xinT=sincRigidTransform(xinT,Tb,0,E.Fof,E.Fob,0);%YB: Indicate here as well that it is an inverse transform
                %DEPHASING - motion dependent
                if isfield(E,'intraVoxDeph') && E.intraVoxDeph && isfield(E,'Db')&&isfield(E,'Dl');xinT=bsxfun(@times,xinT,intravoxDephasing(Tr,dynInd(E.Db.cr,vA(vA<=E.NMs),5),size(xinT), E.Db,E.Dl));            end
                %if isfield(E,'Ds');xinT=bsxfun(@times,xinT,conj(dephaseRotationSusc(dynInd(Tr,E.Ds.d,6),E.Ds.chi,E.Ds.TE)));end   % Susceptibility model with voxel basis
                if isfield(E,'Ds') && isfield(E.Ds,'Dl')% Susceptibility model with voxel basis
                    TrD = dynInd(E.Ds.Tr,vA(vA<=E.NMs),5);
                    xinT=bsxfun(@times,xinT,conj(dephaseRotationTaylor(TrD,E.Ds.f,E.Ds.D,E.Ds.TE)));
                end  
                %if isfield(E,'Dl');xinT=bsxfun(@times,xinT,conj(dephaseRotation(dynInd(Tr,E.Dl.d,6),E.Dl.D,E.Dl.TE)));end
                if isfield(E,'Dl')
                    TrD=dynInd(E.Dl.Tr,vA(vA<=E.NMs),5); 
                    xinT=bsxfun(@times,xinT,conj(dephaseRotationTaylor(TrD, E.Dl.f,E.Dl.D,E.Dl.TE)));
                end
                if isfield(E,'Dc');xinT=bsxfun(@times,xinT,conj(dephaseRotation(dynInd(Tr,E.Dc.d,6),E.Dc.D)));end  
            end
            %DEPHASING - motion independent
            if isfield(E,'Db');xinT=bsxfun(@times,xinT,conj(dephaseBasis( E.Db.B, dynInd(E.Db.cr,vA(vA<=E.NMs),5), dynInd(size(xinT),1:3,2),E.Db.TE)) );end   
            if isfield(E,'B0Exp');xinT=bsxfun(@times,xinT, conj(exp(1i* 2*pi*E.B0Exp.TE * dynInd(E.B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5))));end     
        end
        xouT=dynInd(xouT,vA<=E.NMs,5,xinT);
        
        if notsum ==1 && ~isempty(xouT) % YB: Used in GDsolver.m
            if isempty(xou_notsum); xou_notsum = xouT;
            else xou_notsum=cat(5,xou_notsum,xouT);end 
        end
        if ~isempty(xouT);xouT=sum(xouT,5);end
    end
    if isfield(EH,'Fms') && ~isempty(xouT);xouT=sum(xouT,5);end
            
    %SLICE MASK
    if isfield(E,'Bm') && ~isempty(E.Bm);xouT=sum(xouT,6);end

    if a==E.oS(1) || isempty(xou);xou=xouT;elseif (isfield(E,'Tr') || isfield(EH,'Fms')) && ~isempty(xouT);xou=xou+xouT;elseif ~isempty(xouT);xou=cat(5,xou,xouT);end
end;xT=[];xouT=[];

if isfield(E,'Ds') && isfield(E.Ds,'B0')% Susceptibility model with voxel basis
    xou = bsxfun(@times, xou , conj(exp(1i*2*pi*E.Ds.TE*E.Ds.B0)));
    if ~isempty(xou_notsum);xou_notsum = bsxfun(@times, xou_notsum , conj(exp(1i*2*pi*E.Ds.TE*E.Ds.B0)));end
end

if isfield(EH,'Fb')%Filtering
    x=xou;   
    if isfield(EH.Fb,'Fd')
        for n=1:EH.Fb.ND
            padEl=zeros(1,EH.Fb.ND+1);padEl(n)=1;
            xs=diff(padarray(dynInd(x,n,EH.Fb.ND+1),padEl,0,'pre'),1,n)/EH.Fb.w(n);
            if n==1;xou=xs;else xou=xou+xs;end
        end            
    else
        ND=numDims(x);
        cont=1;
        if length(EH.Fb)>1;extDim=1;else extDim=0;end
        for n=1:length(EH.Fb)        
            mirr=zeros(1,ND-extDim);mirr(n)=EH.ma(n);       
            if EH.ma(n)
                xs=cat(n,dynInd(x,cont,ND),dynInd(x,cont+1,ND));cont=cont+2;
            else
                xs=dynInd(x,cont,ND+1-extDim);cont=cont+1;
            end
            if isa(EH.Fb{n},'gpuArray');xs=gpuArray(xs);end
            xs=filtering(xs,EH.Fb{n},EH.mi(n));
            xs=mirroring(xs,mirr,0);
            if n==1;xou=xs;else xou=xou+xs;end
        end
    end
end

%SLAB EXTRACTION
if isfield(E,'ZSl');xou=extractSlabs(xou,abs(E.ZSl),0,0);end

if isfield(EH,'NXRec');xou=resampling(xou,EH.NXRec);end

if isfield(EH,'Mc');xou=bsxfun(@times,xou,EH.Mc);end%Mask for sensitivity computation
if isfield(EH,'Xb');xou=bsxfun(@times,xou,EH.Xb);end%Body coil or mask
if isfield(EH,'Gb');xou=filtering(xou,EH.Gb,EH.mi);end%Filtering

if isfield(EH,'Ab');xou=EH.Ab*xou;end%Explicit matrix


