
function [Ab_grad, E, x, PT] = gradPTCalibration(y,E,EH, x,C,deb,PT)

%TODO: enable use of RAS motion parameters
%TODO: use EH.We to downweight certain shots in calibration

if nargin<4;C=[];end
if nargin<5 || isempty(deb);deb=2;end

if isfield(E,'Tr')              
    %GENERAL PARAMETERS AND ARRAYS
    gpu=isa(x,'gpuArray');
    T=E.Tr;E=rmfield(E,'Tr'); %YB: need to remove it otherwise the encode will perform the tranfsormation twice
    if isfield(E,'Dl');TD=E.Dl.Tr;end
    if isfield(E,'Ds');TD=E.Ds.Tr;end
    
    if isfield(E,'Db'); Db = E.Db; E = rmfield(E,'Db');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dc'); Dc = E.Dc; E = rmfield(E,'Dc');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dl'); Dl = E.Dl; E = rmfield(E,'Dl');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Ds'); Ds = E.Ds; E = rmfield(E,'Ds');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'B0Exp'); B0Exp = E.B0Exp; E = rmfield(E,'B0Exp');end%YB: need to remove it otherwise the encode will dephase it twice
    if isfield(E,'Dbs'); Dbs = E.Dbs; E = rmfield(E,'Dbs');end%Could include this part in encode.m, but it will be skipped there as it in in the if E.Tr....end
    if isfield(E,'B1m'); B1m = E.B1m; E = rmfield(E,'B1m');end
    if isfield(E,'Shim'); Shim = E.Shim; E = rmfield(E,'Shim');end
    
    NT=size(T);
    NX=size(x);    
    dimG=max(6,length(NT));dimS=find(NT(1:dimG-1)~=1);if isempty(dimS);dimS=dimG-1;end%Dimensions of motion states and parameters of the transform 
    ndS=NT(dimS);ndG=NT(dimG);NTS=length(dimS);
    ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
    ndH=size(ha,2);   
    flagw=zeros(NT(1:dimS(end)));
    En=single(zeros([ndS 1]));EnPrev=En;EnPrevF=EnPrev;EnF=En;   
    dH=single(zeros([ndS ndH]));dG=single(zeros([ndS ndG]));dGEff=dG;    %YB: For every state you have 21 Hessian elements and 6 gradient terms    
    multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)    
    NElY=numel(y);
    if isfield(E,'nF');E.Sf=dynInd(E.Sf,E.nF,3);y=dynInd(y,E.nF,3);end%Extract ROI in the readout direction
    if isfield(E,'Ps') && ~isempty(E.Ps);E.Sf=bsxfun(@times,E.Sf,E.Ps);y=bsxfun(@times,y,E.Ps);end%Preconditioner of the coils, deprecated
    
    %BLOCK SIZES AND CONVERGENCE VALUES FOR THE TRANSFORM
    ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;ET.kG=E.kG;ET.rkG=E.rkG;   
    if isfield(E,'Bm');ET.Bm=E.Bm;ET.bS(1)=1;end%YB: in encode.m = slice mask
    if isfield(E,'Fms');ET.Sf=E.Sf;ET.Fms=E.Fms;ET.bS(1)=NX(3);E.dS(1)=1;end%YB: in encode.m = Fourier encoding multislice
    if isfield(E,'ZSl');ET.ZSl=E.ZSl;E.ZSl=-E.ZSl;end%YB: in encode.m = Slab extraction
    ET.cTFull=E.cT;       

    if NTS==1;NB=1;else NB=NT(dimS(1));end
    for b=1:NB
        if NB==1;ET.cT=ET.cTFull;else ET.cT=dynInd(ET.cTFull,b,dimS(1));end
        if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

        %COMPUTE THE JACOBIAN AND STORE VALUES
        ind2EstOr=find(~ET.cT(:));NEst=length(ind2EstOr); %YB: indeces which state to update (so not converged)
        if isfield(E,'Fs');ry=zeros([size(y,1) 1],'like',real(y));else ry=zeros([NT(6) 1],'like',real(y));end

        if exist('Ds','var') && isfield(Ds,'B0')% Susceptibility model with voxel basis
            x = bsxfun(@times, x , exp(1i*2*pi*Ds.TE*Ds.B0));
        end

        for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2EstOr(vA);
            if isfield(E,'Fms')
               [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
               E.Sf=dynInd(ET.Sf,vA,3);            
            end                               
            if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
            if NB==1;vT=vA;else vT={b,vA};end
            Ti=dynInd(T,vT,dimS);
            [Tf,Tg]=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,1,1,1);%Transform factors 

            xT=x;
            if exist('B0Exp','var');xT=bsxfun(@times,xT, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end  
            if exist('Db','var');xT=bsxfun(@times,xT,dephaseBasis( Db.B, dynInd(Db.cr,vA(vA<=E.NMs),5), size(xT),Db.TE));end
            if exist('Dc','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
            %if exist('Dl','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dl.d,6),Dl.D,Dl.TE));end   % Linear model with voxel model
            if exist('Dl','var')
                TiD = dynInd(TD,vT,dimS);
                xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Dl.f,Dl.D,Dl.TE));
            end
            if exist('Ds','var') && isfield(Ds,'D')
                TiD = dynInd(TD,vT,dimS);
                xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Ds.f,Ds.D,Ds.TE));
            end   % Susceptibility model with voxel basis
            if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xT=bsxfun(@times,xT,intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xT),Db,Dl));end
              
            if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
            if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                       
            [xT,xTG]=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);
            
            if exist('Dbs','var');xT=bsxfun(@times,xT,dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xT),Dbs.TE));end
            if exist('B1m','var');xT=bsxfun(@times,xT,exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xT),[],1)));end   
            if exist('Shim','var');xT=bsxfun(@times,xT,exp(+1i*2*pi*Shim.TE *Shim.B0) );end   

            if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end%YB: Taking out E.nF only after transformation, otherwise transformations not consistent with rest of algorithm
            if isfield(E,'Fs')
                E.vA=vA;
                E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                xT=encode(xT,E)-dynInd(y,indY,1); % YB:This is to create the term w in equation 11 of Lucilio's first paper
            elseif isfield(E,'Bm')
                E.Bm=dynInd(ET.Bm,vA,6);
                xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
            else
                xT=encode(xT,E)-dynInd(y,vT,[5 3]);
            end
            
            if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
            if isfield(E,'Fs')
                ry(indY)=normm(xT,[],3:4)/NElY;
                EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
            else
                ry(vA)=normm(xT)/NElY;
                EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(ry(vA)));
            end
            if isfield(E,'Fd') && ~isempty(E.Fd)%Filtering
                xT=fftGPU(xT,3);
                xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                EnPrev=dynInd(EnPrev,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
            else
                EnPrev=dynInd(EnPrev,vT,1:NTS,dynInd(EnPrevF,vT,1:NTS));
            end

            xTG=sincRigidTransformGradient(xTG,Tf,Tg,E.Fof,E.Fob); % This creates outer 2 operations of w_l in Lucilios 1st paper (still needs encoding!)

            for g=1:ndG %YB: Building Gradient (with the use of the Jacobian)
                %Apply phase in scanner reference to Gradient as well
                if exist('Dbs','var');xTG{g}=bsxfun(@times,xTG{g},dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xTG{g}),Dbs.TE));end
                if exist('B1m','var');xTG{g}=bsxfun(@times,xTG{g},exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xTG{g}),[],1)));end   
                if exist('Shim','var');xTG{g}=bsxfun(@times,xTG{g},exp(+1i*2*pi*Shim.TE *Shim.B0) );end 
                
                if isfield(E,'nF');xTG{g}=dynInd(xTG{g},E.nF,3);end
                xTG{g}=encode(xTG{g},E); % YB: This is the encoding step still needed to become term w_l in 11th equation in 1st paper
                if isfield(E,'Pf') && ~isempty(E.Pf);xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Pf,indY,1));end  %Precondition using all the residuals, deprecated
                if isfield(E,'Fd') && ~isempty(E.Fd)
                    xTG{g}=fftGPU(xTG{g},3);
                    xTG{g}=bsxfun(@times,xTG{g},dynInd(E.Fd,indY,1));
                end
%                 %%% YB extra dephasing term to gradient
%                 if exist('Dl','var') && size(Dl.D,6)==2
%                     if ismember(g, Dl.d)
%                         E.Tr = T; E.Dl = Dl; if exist('Db','var');E.Db=Db;end % Must be T since Ti will be extracted via bS in encode
%                         %xD = dynInd(Dl.D,g==Dl.d,6).* x;%dynInd works with logical indices
%                         xD = sum( bsxfun(@times, x, dynInd(Dl.fDeriv(E.Tr), g==[4:6], 7) ) , 6); 
%                         if isfield(E,'nF')
%                             xD=dynInd(xD,E.nF,3);E.Dl.D=dynInd(E.Dl.D,E.nF,3);
%                             if exist('Db','var')
%                                 B = E.Db.B; B = reshape(B, [NX, size(B,2)]);B=dynInd(B,E.nF,3);
%                                 sizeB=size(B);B = reshape(B,[prod(sizeB(1:3)),size(B,4)]);E.Db.B = B;B=[];
%                             end
%                         end %D field must also be changed temporarily
%                         %xTG_Linear = 2*pi*1i*Dl.TE *(180/pi) * encode(xD,E);%180/pi because dephaseRotation in degrees
%                         xTG_Linear = 2*pi*1i*Dl.TE * encode(xD,E);%NO 180/pi because dephaseRotation in radians
%                         xTG{g} = xTG{g} + xTG_Linear;%Modification to equation 11 in AlignedReconstruction
%                         E = rmfield(E,{'Tr','Dl'});
%                         if isfield(E,'Db'); E = rmfield(E,'Db');end
%                     end
%                 end
                
                if isfield(E,'Fs')
                    dG(vA,g)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xTG{g}.*conj(xT)),3:4)),vA,1));
                elseif isfield(E,'Bm')
                    dG(vA,g)=gather(multDimSum(real(xTG{g}.*conj(xT)),1:6));
                else
                    dG=dynInd(dG,[vT g],1:NTS+1,gather(permute(multDimSum(real(xTG{g}.*conj(xT)),[1:2 4]),[5 3 1 2 4])));
                end
            end

        end;xTG=[];
    end

    %%% TO DO: convert gradients to RAS frame if needed
    
    %%% CONVERT MOTION GRADIENTS TO CALIBRATION GRADIENTS   - this is the new part vs
    %dG=restrictTransform(dG);%Rotation parameters between -pi and pi
    Ab_grad = zerosL(PT.Ab);
    dimStates=1;
    dimMotParam=2;
    for ii=1:size(Ab_grad,1)%Motion params
        for jj=1:size(Ab_grad,2)%NCha for each PT signal
           %%% Chain rule
           gradTemp = dynInd(dG, ii, dimMotParam);
           if jj==size(Ab_grad,2) && jj>size(PT.pTimeResGr,1)
               ptSignal=1;%The offset
           else
               ptSignal = dynInd(PT.pTimeResGr,jj,1).';%Extract certain coil
           end
           gradTemp = gradTemp .* ptSignal;
           
           Ab_grad(ii,jj) = gather(multDimMea(gradTemp,dimStates)) ; %Sum over states
        end
    end
    
    
    %%% REST
    T=constrain(T,C);
    E.Tr=T;E.cT=ET.cTFull;E.bS=ET.bS;E.dS=ET.dS;E.oS=ET.oS;E.kG=ET.kG;E.rkG=ET.rkG;
    if exist('Dl','var');Dl.Tr = TD;end%Updated one
    if exist('Ds','var') && isfield(Ds, 'D');Ds.Tr = TD;end%Updated one
    if exist('B0Exp','var'); E.B0Exp = B0Exp;end%YB add again
    if exist('Db','var'); E.Db = Db;end%YB add again
    if exist('Dc','var'); E.Dc = Dc;end%YB add again
    if exist('Dl','var'); E.Dl = Dl;end%YB add again
    if exist('Ds','var'); E.Ds = Ds;end%YB add again
    if exist('Dbs','var'); E.Dbs = Dbs;end%YB add again
    if exist('B1m','var'); E.B1m = B1m;end%YB add again
    if exist('Shim','var'); E.Shim = Shim;end%YB add again

    if isfield(E,'Fms');E.Sf=ET.Sf;E.Fms=ET.Fms;end
    if isfield(E,'ZSl');E.ZSl=ET.ZSl;end
    if isfield(E,'Bm');E.Bm=ET.Bm;end
    if isfield(E,'vA');E=rmfield(E,'vA');end
    if isfield(E,'Fd');E=rmfield(E,'Fd');end   
    Enmin=min(EnPrev,En);
    if ~isfield(E,'En');E.En=Enmin;else E.En(ind2EstOr)=Enmin(ind2EstOr);end%To keep a record of the achieved energy
end
