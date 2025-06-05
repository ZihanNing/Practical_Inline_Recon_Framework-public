function [svr,conv]=svrSolveTVolu(svr)

%SVRSOLVET   Solve for T on a per-volume level
%   SVR=SVRSOLVET(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

NT=size(svr.TV);ndT=ndims(svr.TV);    

mirr=[2 2 2];
a=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;
   1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];
NHe=size(a,2);
dHe=single(zeros([NHe NT(ndT-1)]));
dH=single(zeros(NT([ndT ndT-1])));dHEff=dH;
E=single(zeros(NT(1:ndT-1)));
if ~isfield(svr,'convT');convT=single(false(NT(1:ndT-1)));else convT=svr.convT;end
if isfield(svr,'EnViews');convT(~svr.EnViews)=1;end

multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)
     
%winic=1e-2;
winic=svr.ParSVR.winic;
if ~isfield(svr,'w');svr.w=winic*ones(NT(1:ndT-1));end
svr.w=min(svr.w,winic);
flagw=zeros(NT(1:ndT-1));
perm=1:ndT;perm([1 2 ndT-1 ndT])=[ndT-1 ndT 1 2];

Eprev=E;
dHe(:)=0;dH(:)=0;

%TRANSFORMING
[et,etg]=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,svr.TV,1,1,[],1,svr.cGridV);
[xT,xB]=sincRigidTransform(svr.x,et,1,svr.FTV,svr.FTHV);xT=real(xT);
MT=real(sincRigidTransform(svr.M,et,1,svr.FTV,svr.FTHV));
G=sincRigidTransformGradient(xB,et,etg,svr.FTV,svr.FTHV);

GT=cell(1,NT(ndT));GTH=cell(1,NT(ndT));
vT=find(~convT)'; 
for v=vT     
    [svr.yy{v},svr.MM{v}]=computeResiduals; 
    Eprev(v)=gather(multDimSum(real(svr.yy{v}.*conj(svr.yy{v})),[1:ndT-2 5]));

    for m=1:NT(ndT)
        GT{m}=encode(real(G{m}));
        GT{m}=filterMask(GT{m});       
        GTH{m}=conj(GT{m});
    end

    for m=1:NHe;dHe(m,v)=gather(multDimSum(real(GT{a(1,m)}.*GTH{a(2,m)}),[1:ndT-2 5]));end   
    for m=1:NT(ndT);dH(m,v)=gather(multDimSum(real(GTH{m}.*svr.yy{v}),[1:ndT-2 5]));end
end

E=Eprev;
MHe=single(eye(NT(ndT)));
flagw(:)=0;    
fina=0;
while fina==0
    vT=find(~convT)';
    for s=vT
        for k=1:NHe
            if a(1,k)==a(2,k)              
                MHe(a(1,k),a(2,k))=(1+svr.w(s))*dHe(k,s)+1e-9;%1e-9 serves to stabilize
            else
                MHe(a(1,k),a(2,k))=dHe(k,s);MHe(a(2,k),a(1,k))=dHe(k,s);
            end              
        end   
        dHEff(:,s)=-winic*single(double(MHe)\double(dH(:,s)))/svr.w(s);
    end     
    dHEff(:,svr.w>1e10 | convT)=0;
    permH=1:ndT;permH(ndT-1:ndT)=[2 1];permH(1:ndT-2)=3:ndT;
    Tupr=permute(dHEff,permH);
    Tup=svr.TV+Tupr;
    Tup=restrictTransform(Tup);    
    et=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,Tup,1,0,[],1,svr.cGridV);                   
    xT=real(sincRigidTransform(svr.x,et,1,svr.FTV,svr.FTHV));
    MT=real(sincRigidTransform(svr.M,et,1,svr.FTV,svr.FTHV));

    for v=vT
        [svr.yy{v},svr.MM{v}]=computeResiduals;
        E(v)=gather(multDimSum(real(svr.yy{v}.*conj(svr.yy{v})),[1:ndT-2 5]));
    end
    E(svr.w>1e10 & ~convT)=Eprev(svr.w>1e10 & ~convT);
    flagw(E<=Eprev)=2;               
    %fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(Eprev),sum(E));

    if any(flagw==1 | flagw==0)  
        svr.w(E>Eprev & ~convT)=svr.w(E>Eprev & ~convT)*multB;
        %svr.w(E>Eprev)=svr.w(E>Eprev)*multB;
    else   
        svr.w(~convT)=svr.w(~convT)/multA;
        %svr.w=svr.w/multA;
        svr.w(svr.w<1e-8)=multA*svr.w(svr.w<1e-8);%To avoid numeric instabilities 
        %Tupr=bsxfun(@minus,Tupr,mean(Tupr,ndT-1));%This gives problems with regularization, so better to disable
        svr.TV=svr.TV+Tupr;
        svr.TV=restrictTransform(svr.TV);
        fina=2;
        traMax=abs(permute(dynInd(Tupr,1:3,ndT),perm));
        rotMax=convertRotation(abs(permute(dynInd(Tupr,4:6,ndT),perm)),'rad','deg');
        fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(Eprev),sum(E));
%         fprintf('Maximum change in translation (vox): ');fprintf('%0.3f ',max(traMax,[],1));
%         fprintf('/ Maximum change in rotation (deg): ');fprintf('%0.3f ',max(rotMax,[],1));fprintf('\n');
        fprintf('Maximum change in translation (vox): %s/ rotation (deg): %s\n',sprintf('%0.3f ',max(traMax,[],1)),sprintf('%0.3f ',max(traMax,[],1)));
        traLim=0.16;rotLim=0.08;
        traLim=traLim*svr.ParSVR.convL;rotLim=rotLim*svr.ParSVR.convL;
        if max(traMax(:))>traLim || max(rotMax(:))>rotLim;fina=1;end
        convT(max(traMax,[],2)<traLim & max(rotMax,[],2)<rotLim)=1;
        conv=all(convT);            
        fprintf('Not converged motion states: %s of %s\n',num2str(NT(ndT-1)-sum(single(convT))),num2str(NT(ndT-1)));
        if conv;svr=rmfield(svr,'w');end
    end 
end        

if ~isfield(svr,'TVHist');svr.TVHist=svr.TV;else svr.TVHist=cat(1,svr.TVHist,svr.TV);end
svr.convT=convT;
if conv;svr=rmfield(svr,'convT');end

function [x,M]=computeResiduals
    x=encode(xT);
    svr.MM{v}=encode(MT);        
    x=bsxfun(@minus,x,svr.y{v});
    x=filterMask(x);
    M=svr.MM{v};
end

function x=encode(x)
    x=x(:,:,:,v);       
    x=resampling(x,svr.NN,2);
    xa=zeros(svr.NmReal(1,:,v),'like',x);       
    xa=dynInd(xa,svr.vAcq{v},1:3,dynInd(x,svr.vRec{v},1:3));        
    if ~svr.ParSVR.quickRecon;x=resampling(xa,svr.NYReal(1,:,v),[],mirr);%%%%THIS TAKES TIME, ACCELERATED BELOW
    elseif svr.ParSVR.quickRecon==1;x=resamplingQuick(xa,svr.NYReal(1,:,v));
    else x=resamplingQuick(xa,svr.NYReal(1,:,v),'nearest');
    end
    if svr.ParSVR.quickRecon~=3;x=filtering(x,svr.SlPr{v},1);end
end

function x=filterMask(x)
    if svr.ParSVR.fracOrd~=0;x=filtering(x,svr.G{v},[],[],1,svr.FY{v},svr.FYH{v});end%%%%%NOTE WE'VE INTRODUCED THE QUICK FLAG TO ACCELERATE IT
    x=bsxfun(@times,x,max(svr.MM{v},1e-4));
    %x=bsxfun(@times,x,svr.MM{v});%IT SHOULD BE SQRT WT, THIS IS A DIRTY TRICK TO MAKE IT WORK BETTER... ONCE TRACKING IMPLEMENTED PROBABLY NOT NECESSARY ANYMORE
end

end