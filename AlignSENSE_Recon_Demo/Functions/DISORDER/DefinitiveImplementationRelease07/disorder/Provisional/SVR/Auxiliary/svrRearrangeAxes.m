function svr=svrRearrangeAxes(svr)

%SVRREARRANGEAXES   Sets up the axes for SVR: 
%Rearrange the images so that the arrays are as close as possible to the
%reference image, store new slice and in plane dimensions
%   SVR=SVRREARRANGEAXES(SVR)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   * SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

svr.NV=length(svr.y);
svr.id=zeros(svr.NV,3);%Dimensions after permuting the axes
for v=1:svr.NV
    MT=svr.MT{v}(1:3,1:3);
    MTn=bsxfun(@rdivide,MT,sqrt(sum(MT.^2,1)));%Orientation vector along different axis
    [~,svr.id(v,:)]=max(abs(MTn),[],1);%Closest material dimension for each spatial dimension
    ord=sign(indDim(MTn,svr.id(v,:),1));%Direction of each spatial dimension
    for n=1:3%We flip the dimensions
        if ord(n)<0
            svr.y{v}=flip(svr.y{v},n);
            svr.yOr{v}=flip(svr.yOr{v},n);
            svr.El{v}=flip(svr.El{v},n);
            svr.MPack{v}=flip(svr.MPack{v},n);
            for p=1:svr.P(v)
                svr.MSlices{v}{p}=flip(svr.MSlices{v}{p},n);
                svr.MSlicesP{v}{p}=flip(svr.MSlicesP{v}{p},n);                
            end
            Nf=zeros(4,1);Nf(n)=size(svr.y{v},n)-1;Nf(4)=1;
            svr.MT{v}(:,4)=svr.MT{v}*Nf;
            svr.MT{v}(:,n)=-svr.MT{v}(:,n);
            svr.Elp{v}(n)=Nf(n)+2-svr.Elp{v}(n);
        end
    end%We permute the dimensions
    svr.y{v}=ipermute(svr.y{v},svr.id(v,:));
    svr.yOr{v}=ipermute(svr.yOr{v},svr.id(v,:));
    svr.El{v}=ipermute(svr.El{v},svr.id(v,:));
    svr.MPack{v}=ipermute(svr.MPack{v},[svr.id(v,:) 4 5]);
    for p=1:svr.P(v)
        svr.MSlices{v}{p}=ipermute(svr.MSlices{v}{p},[svr.id(v,:) 4 5 6]);
        svr.MSlicesP{v}{p}=ipermute(svr.MSlicesP{v}{p},[svr.id(v,:) 4 5 6]);
    end
    svr.MS{v}(svr.id(v,:))=svr.MS{v};
    svr.MT{v}(1:3,svr.id(v,:))=svr.MT{v}(1:3,1:3);
    svr.Elp{v}([svr.id(v,:) 3+svr.id(v,:)])=svr.Elp{v};
    svr.NY{v}=size(svr.y{v})';%Sizes
end
%svrWriteData(svr,'Pe',svr.y,1,[],'',0);
svr.MS=cat(3,svr.MS{:});svr.MS=permute(svr.MS,[2 1 3]);
svr.MT=cat(3,svr.MT{:});
svr.NY=cat(3,svr.NY{:});
svr.Elp=cat(3,svr.Elp{:});