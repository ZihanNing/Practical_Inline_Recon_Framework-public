function svr=svrDecode(svr,sep,fil,arr)

%SVRDECODE   Decoding operator for SVR
%   SVR=SVRDECODE(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<2 || isempty(sep);sep=0;end
if nargin<3 || isempty(fil);fil=1;end
if nargin<4 || isempty(arr);arr=0;end

mirr=[2 2 2];

y=zeros([svr.NNEncode svr.NV],'like',svr.x);
%RECONSTRUCTED DATA
for v=1:svr.NV     
    if ~isfield(svr,'EnViews') || svr.EnViews(v)

        %ARRANGING TO REMOVE DROPOUT SLICES 
        if arr
            Res=4;
            x=permute(svr.yy{v},svr.id(v,:));
            NXA=size(x);NXA(end+1:3)=1;
            x=resampling(x,[NXA(1:2) Res*ceil(NXA(3)/Res)],2);
            xM=multDimMea(abs(x),1:2);
            NX3=size(x,3);
            xM=reshape(xM,[Res NX3/Res]);
            [~,iMax]=max(xM,[],1);
            nMax=(1:Res:NX3)+iMax-1;        
            x=x(:,:,nMax);
            x=repmat(x,[1 Res 1]);
            x=reshape(x,[NXA(1:2) NX3]);
            x=resampling(x,NXA,2);
            x=ipermute(x,svr.id(v,:));          
        else
            x=svr.yy{v};
        end
        
        %SLICE WEIGHTS
        x=bsxfun(@times,x,svr.W{v});

        %APODIZATION
        if svr.ParSVR.UseApo>0;x=bsxfun(@times,x,svr.H{v});end

        %PER-PACK MASKING
        if any(svr.TP{v}(:)~=0);x=bsxfun(@times,x,svr.MPack{v});end

        for p=1:size(x,5)
            xp=x(:,:,:,:,p);

            %PER-SLICE MASKING
            if any(svr.TE{v}{p}(:)~=0);xp=bsxfun(@times,xp,svr.MSlices{v}{p});end

            %SLICE PROFILE        
            if fil && svr.ParSVR.SlBefore==0;xp=filtering(xp,svr.SlPr{v},1);end

            %RESAMPLING
            xp=resampling(xp,svr.NmReal(1,:,v),[],mirr);

            %SLICE PROFILE
            if fil && svr.ParSVR.SlBefore==1;xp=filtering(xp,svr.SlPr{v},1);end

            %FILL ARRAY
            NNXX=size(xp);NNXX(end+1:4)=1;    
            z=zeros([svr.NN NNXX(4:end)],'like',svr.x);
            z=dynInd(z,svr.vRec{v},1:3,dynInd(xp,svr.vAcq{v},1:3));
            z=resampling(z,svr.NNEncode,2);

            if any(svr.TE{v}{p}(:)~=0)
                et=precomputeFactorsSincRigidTransform(svr.kGrid,svr.rkGrid,svr.TE{v}{p},0,0,[],1,svr.cGrid);
                z=sum(real(sincRigidTransform(z,et,1,svr.FT,svr.FTH)),6);%THIS LOOKS LIKE A BUG, AS IT SHOULD BE DI=0 PROBABLY!!!
            end
            if p==1;zo=z;else zo=cat(5,zo,z);end
        end

        if any(svr.TP{v}(:)~=0)
            %PER-PACK TRANSFORM
            et=precomputeFactorsSincRigidTransform(svr.kGrid,svr.rkGrid,svr.TP{v},0,0,[],1,svr.cGrid);
            y(:,:,:,v)=sum(real(sincRigidTransform(zo,et,0,svr.FT,svr.FTH,0)),5);
        else
            y(:,:,:,v)=zo;                    
        end
    end
end
if isfield(svr,'EnViews') && sep==0;y=dynInd(y,svr.EnViews,4);end
y=resampling(y,svr.NN,2);

if sep==2
    svr.xx=y;
else
    %TRANSFORMING
    if isfield(svr,'EnViews') && sep==0;et=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,dynInd(svr.TV,svr.EnViews,4),0,0,[],1,svr.cGridV);
    else et=precomputeFactorsSincRigidTransform(svr.kGridV,svr.rkGridV,svr.TV,0,0,[],1,svr.cGridV);
    end
    y=real(sincRigidTransform(y,et,0,svr.FTV,svr.FTHV,0));
    if sep==0;svr.xx=sum(y,4);elseif sep==1;svr.xx=y;end
end
svr.xx=real(svr.xx);
