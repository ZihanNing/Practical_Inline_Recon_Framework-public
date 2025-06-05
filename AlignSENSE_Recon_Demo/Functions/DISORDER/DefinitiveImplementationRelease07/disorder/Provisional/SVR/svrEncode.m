function svr=svrEncode(svr,we)

%SVRENCODE   Encoding operator for SVR
%   SVR=SVRENCODE(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   ** SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2 || isempty(we);we=1;end

%TRANSFORMING
if isfield(svr,'EnViews');et=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,dynInd(svr.TV,svr.EnViews,4),1,0,[],1,svr.cGridV);
else et=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,svr.TV,1,0,[],1,svr.cGridV);
end
y=real(sincRigidTransform(svr.xx,et,1,svr.FTV,svr.FTHV));
y=resampling(y,svr.NNEncode,2);

mirr=[2 2 2];

indV=1:svr.NV;
if isfield(svr,'EnViews');indV=indV(svr.EnViews);end
c=1;
for v=1:svr.NV
    if ~isfield(svr,'EnViews') || svr.EnViews(v)==1    
        x=y(:,:,:,c);
    
        %PER-PACK TRANSFORM
        if any(svr.TP{indV(c)}(:)~=0)
           et=precomputeFactorsSincRigidTransform(svr.kGrid,svr.rkGrid,svr.TP{indV(c)},1,0,[],1,svr.cGrid);
           x=real(sincRigidTransform(x,et,1,svr.FT,svr.FTH));
        end

        for p=1:size(x,5)        
            xp=x(:,:,:,:,p);
            %PER-SLICE TRANSFORM
            if any(svr.TE{indV(c)}{p}(:)~=0)
                et=precomputeFactorsSincRigidTransform(svr.kGrid,svr.rkGrid,svr.TE{indV(c)}{p},1,0,[],1,svr.cGrid);
                xp=real(sincRigidTransform(xp,et,1,svr.FT,svr.FTH));
            end   
            %EXTRACT ARRAYS
            xp=resampling(xp,svr.NN,2);  
            NNXX=size(xp);NNXX(end+1:4)=1;
            xa=zeros([svr.NmReal(1,:,indV(c)) NNXX(4:end)],'like',xp);        
            xa=dynInd(xa,svr.vAcq{indV(c)},1:3,dynInd(xp,svr.vRec{indV(c)},1:3));     

            %SLICE PROFILE
            if svr.ParSVR.SlBefore==1;xa=filtering(xa,svr.SlPr{indV(c)},1);end

            %RESAMPLING
            xa=resampling(xa,svr.NYReal(1,:,indV(c)),[],mirr); 

            %SLICE PROFILE
            if svr.ParSVR.SlBefore==0;xa=filtering(xa,svr.SlPr{indV(c)},1);end

            %PER-SLICE MASKING
            if any(svr.TE{indV(c)}{p}(:)~=0);xa=sum(bsxfun(@times,xa,svr.MSlices{indV(c)}{p}),6);end

            if p==1;xo=xa;else xo=cat(5,xo,xa);end
        end
        x=xo;

        %PER-PACK MASKING
        if any(svr.TP{indV(c)}(:)~=0);x=sum(bsxfun(@times,x,svr.MPack{indV(c)}),5);end

        %USE THE WEIGHTS FOR RECONSTRUCTION
        %if we~=0;x=bsxfun(@times,x,svr.W{indV(c)});end    

        %ASSIGNMENT
        svr.yy{indV(c)}=x;
        c=c+1;
    end
end



