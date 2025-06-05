addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06'));
addpath(genpath('/home/lcg13/Work/aloha_public/aloha_public/kt/CPU'));
addpath(genpath('/home/lcg13/Work/aloha_public/aloha_public/kt/CPU/../_data'));
addpath(genpath('/home/lcg13/Work/aloha_public/Low-Rank-Matrix-Completion-master'));

tic
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/DataMBAnthony1');
%load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/DataMBAnthony2');
toc
tic
%x=permute(x,[1 2 3 4 6 5]);
xy=bsxfun(@times,S,x);
xy=aplGPU(F,xy,2);
xy=bsxfun(@times,xy,B);
NY=size(y);NY(end+1:6)=1;
NX=size(xy);NX(end+1:6)=1;
xy=resSub(sum(resSub(xy,3,[NY(3) NX(3)/NY(3)]),4),3:4);
y=permute(y,[5 6 1 2 3 4]);
y=diagm(y);
y=permute(y,[3 4 5 6 1 2]);
xy=permute(xy,[6 5 2 4 1 3]);
y=permute(y,[6 5 2 4 1 3]);

xhat=xy;xhat(:)=0;
for n=25:32%NX(1)-16
    fprintf('Readout %d of %d\n',n,NX(1));
    for m=2;%1:NX(3)
        xn=dynInd(xy,[n m],5:6);
        yn=dynInd(y,[n m],5:6);
        NXX=size(xn);
        %KZ-TIME-KY-COILS
        xn=reshape(xn,[prod(NXX(1:3)) 1 NXX(4)]);
        yn=reshape(yn,[prod(NXX(1:3)) 1 NXX(4)]);
        NXXX=size(xn);
      
        sRank=33;
        xn=patch2hank_complex_single(xn,NXXX(1),NXXX(2),NXXX(3),sRank,1);
        yn=patch2hank_complex_single(yn,NXXX(1),NXXX(2),NXXX(3),sRank,1);
        yn(yn==0)=nan;
        for l=1:1
            [U,D,V]=svd(xn,'econ');
            U=U*D;                
            opts.init=1;        
            opts.X=U(:,1:min(sRank,size(U,2)));
            opts.Y=V(:,1:min(sRank,size(U,2)))';
            %visReconstruction(yn)
            tic
            [U,V,~]=lmafit_mc_adp_v2_single(size(xn,1),size(xn,2),33,find(~isnan(yn)),yn(~isnan(yn)),opts);
            xn=U*V;
            toc
            tic
            xn=completion(yn);
            toc
            xn(~isnan(yn))=yn(~isnan(yn));
        end
        %ynhat(yn~=0)=yn(yn~=0);
        
        %ynhat=U*V;
        %visReconstruction(yn)
        xn=hank2patch_complex_single(xn,NXXX(1),NXXX(2),NXXX(3),sRank,1);
        xn=reshape(xn,NXX);
        xhat=dynInd(xhat,[n m],5:6,xn);
    end
end
xhat=permute(xhat,[5 3 6 4 2 1]);
xhat=repmat(xhat,[1 1 NX(3)/NY(3)]);
xhat=sum(bsxfun(@times,xhat,conj(B)),6);
xhat=aplGPU(F',xhat,2);
xhat=bsxfun(@rdivide,sum(bsxfun(@times,conj(S),xhat),4),(sum(abs(S).^2,4)+1e-9));
visReconstruction(x)
visReconstruction(xhat)

return


xy=bsxfun(@times,xy,B);

%xy=resSub(sum(resSub(xy,3,[NY(3) NX(3)/NY(3)]),4),3:4);
%xy=repmat(xy,[1 1 NX(3)/NY(3)]);
y=repmat(y,[1 1 NX(3)/NY(3)]);
y=bsxfun(@times,y,conj(B));
visReconstruction(xy)
visReconstruction(y)


size(xy)
size(y)
toc