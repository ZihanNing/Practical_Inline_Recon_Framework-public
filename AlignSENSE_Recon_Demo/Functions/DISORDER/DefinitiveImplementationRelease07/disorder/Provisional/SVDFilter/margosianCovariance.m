function cov=margosianCovariance(Enc,gS,filt,gpu)

% MARGOSIANCOVARIANCE calculates the covariance corresponding to a given Margosian filter
%   [EIGPV,WEIPV]=MARGOSIANCOVARIANCE(ENC,GS,{FILT},{GPU})
%   ENC is the encoding structure. Following fields are used: KRANGE, the 
%   sampled spectral areas; ACQSIZE, the equivalent size of full spectrum.
%   {GS} is the spatial noise amplification. It defaults to 1
%   {FILT} is the type of filter for covariance calculation. Possibilities are
%   'ramp' for a ramp filter, 'zefi' for a zero filled filter or 'iden' to
%   disable computation based on the encoding structure and using only the
%   g-factors. It defaults to 'ramp'
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   * COV is a structure containing the eigenvalues and weights of the 
%   population spectral distribution. It should contain the following 
%   fields:
%       - COV.DIMPE, the PE dimensions for different filters
%       - COV.LAMBDA, the covariance matrix (first two dimensions) over the
%       voxels not belonging to PE dimensions (third dimension) and 
%       different filters (fourth dimension)
%       - COV.ENC, the encoding structure used to generate the population
%       spectral distribution
%

if nargin<2 || isempty(gS);gS=1;end
if nargin<3 || isempty(filt);filt='ramp';end
if nargin<4 || isempty(gpu);gpu=single(gpuDeviceCount && ~blockGPU);end

if gpu;gS=gpuArray(gS);end
cov.Enc=Enc;
NDims=length(Enc.kRange);
Enc.N;Enc.N(end+1:NDims)=1;Enc.N(4:end)=[];

if Enc.FOVSize(3)==1 && Enc.N(3)>1;Enc.FOVSize(3)=Enc.N(3);end
kGridA=generateGrid(Enc.FOVSize,gpu,Enc.FOVSize,ceil((Enc.FOVSize+1)/2));
kGridB=generateGrid(Enc.N,gpu,Enc.N,ceil((Enc.N+1)/2));

indMA=cell(NDims,3);indMB=cell(NDims,3);
Ek=size(unique(cat(2,Enc.kRange{1},Enc.kRange{2}),'rows'),1);
%This is to modify for sequences like ZEBRA, where the first echo is half
%scan and the remaining are fully sampled
flagZebra=(size(Enc.MPSOffcentres,3)==1);%Only activated if one PE, multi-echo + alternating PE not contemplated!
if ~flagZebra
    E=Ek;
else
    if size(unique(Enc.PartFourier,'rows'),1)==1;E=max(Ek,1);else E=max(Ek,size(Enc.PartFourier,1));end
    if E~=Ek
        for e=Ek+1:E
            for m=1:3
                NKP=round((diff(sort(Enc.kRange{m}(1,:)))+1)/Enc.PartFourier(1,m));
                assert(NKP==Enc.AcqSize(1,m),'Dimensions of spectrum without PF (%d) and acquisition dimensions (%d)  do not match',NKP,Enc.AcqSize(1,m));
                if Enc.kRange{m}(1,1)==floor(Enc.AcqSize(1,m)/2);diF=2;else diF=1;end
                NDiscard=round(NKP*(1-Enc.PartFourier(e,m)));
                Enc.kRange{m}(e,diF)=(-1)^(single(diF==1))*floor((Enc.AcqSize(1,m)-single(diF==2))/2)+((-1)^(single(diF==2)))*NDiscard;
                Enc.kRange{m}(e,3-diF)=(-1)^(single(diF==2))*floor((Enc.AcqSize(1,m)-single(diF==1))/2);
            end
            if e~=1;Enc.AcqSize(e,:)=Enc.AcqSize(e-1,:);end
        end
    end
end

invert=1e-3;
gS=gS.^2;
if size(gS,4)~=E;gS=repmat(gS,[1 1 1 E]);end
cov.dimPE=single(zeros([1 E]));%PE dimension (1 and 2)
cov.dimNoPE=single(zeros([2 E]));%PE dimension (1 and 2)

if strcmp(filt,'iden')
    cov.dimPE(:)=2;
    cov.dimNoPE=repmat([1;3],[1 E]);
    perm=[2 1 3 4];
    gS=permute(gS,perm);
else
    gSaux=[];
    for e=1:E
        noPE=0;
        for m=1:NDims
            if ~flagZebra;NR=diff(sort(Enc.kRange{m}(e,:)))+1;
            else NR=diff(sort(Enc.kRange{m}(1,:)))+1;%Changed to contemplate ZEBRA
            end
            if NR~=1 && NR~=Enc.AcqSize(e,m)
                assert(noPE==0,'More than one dimension for covariance calculation: probably too memory demanding');
                cov.dimPE(e)=m;
                perm=1:3;perm(m)=1;perm(1)=m;
                cov.dimNoPE(:,e)=perm(2:3);
                if isempty(gSaux);gSaux=permute(dynInd(gS,e,4),perm);else gSaux=cat(4,gSaux,permute(dynInd(gS,e,4),perm));end           
            end
        end
    end
    if any(cov.dimPE~=0);gS=gSaux;gSaux=[];end
end
if any(cov.dimPE~=0)
    N=size(gS);N(end+1:4)=1;
    gS=reshape(gS,[1 N(1) prod(N(2:3)) N(4)]);
    cov.Lambda=diagm(complex(gS)); 
else
    cov.Lambda=complex(gS);
end
if ~strcmp(filt,'iden')
    for e=1:E
        if cov.dimPE(e)~=0
            m=cov.dimPE(e);
            NR=diff(sort(Enc.kRange{m}(e,:)))+1;

            pfFactor=NR/Enc.AcqSize(e,m);
            cutOffInd=round(Enc.FOVSize(m)*(1-pfFactor)+1);
            cutOff=gather(abs(kGridA{m}(cutOffInd)-0.5));

            indMA{m}{1}=(kGridA{m}(:)<=-cutOff);indMA{m}{2}=(kGridA{m}(:)>-cutOff & kGridA{m}(:)<=cutOff);indMA{m}{3}=(kGridA{m}(:)>cutOff);
            regSamp=Enc.N(m)-Enc.FOVSize(m);
            assert(regSamp>=0,'Reconstruction spectrum is smaller than acquired spectrum');
            fillRegSamp=[ceil(regSamp/2) floor(regSamp/2)];
            if mod(Enc.FOVSize(m),2)==0;flip(regSamp);end
            indMB{m}{1}=vertcat(true(fillRegSamp(1),1),indMA{m}{1},true(fillRegSamp(2),1));
            indMB{m}{2}=vertcat(false(fillRegSamp(1),1),indMA{m}{2},false(fillRegSamp(2),1));
            indMB{m}{3}=vertcat(false(fillRegSamp(1),1),indMA{m}{3},false(fillRegSamp(2),1));

            kyRM=cell(1,3);for n=1:3;kyRM{n}=gather(kGridB{m}(indMB{m}{n}));end
            sz=ones(1,NDims);sz(m)=Enc.N(m);sz(end+1:2)=1;        
            H=single(zeros(sz));
            if strcmp(filt,'ramp');H(indMB{m}{1})=invert;H(indMB{m}{2})=1+(1-invert)*kyRM{2}/cutOff;H(indMB{m}{3})=2-invert;
            elseif strcmp(filt,'zefi');H(indMB{m}{1})=invert;H(indMB{m}{2})=1;H(indMB{m}{3})=1;end
            if gpu;H=gpuArray(H);end
            if abs(Enc.kRange{m}(e,1))>abs(Enc.kRange{m}(e,2))
                H=flip(H,m);
                for n=1:3;indMB{m}{n}=flip(indMB{m}{n});end
            end    
            H=ifftshift(H,m);
            for n=1:3;indMB{m}{n}=ifftshift(indMB{m}{n});end%Not necessary now but it may be if we connect spectral and spatial covariance estimation
            H=diag(complex(H));
            [F,FH]=buildStandardDFTM(N(1),1,gpu);
            Fl=FH{1}*H*F{1};
            Fr=FH{1}*conj(H)*F{1};
            C=dynInd(cov.Lambda,e,4);
            if gpu
                C=pagefun(@mtimes,Fl,C);
                C=pagefun(@mtimes,C,Fr);
            else
                for n=1:size(C,3)
                    C=dynInd(C,n,3,Fl*dynInd(C,n,3));
                    C=dynInd(C,n,3,dynInd(C,n,3)*Fr);
                end
            end
            cov.Lambda=dynInd(cov.Lambda,e,4,C);
        end
    end
end
cov.Lambda=gather(cov.Lambda);

%To visualize the covariance
%for n=1:3
%    visReconstruction([real(cov.Lambda(:,:,300,n)) imag(cov.Lambda(:,:,300,n))],[],[],0);
%end
%pause
