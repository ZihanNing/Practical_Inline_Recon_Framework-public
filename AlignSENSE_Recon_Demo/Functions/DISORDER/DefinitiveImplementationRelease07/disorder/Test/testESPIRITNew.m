tic
%if ~exist('rec','var')
load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/Hankel.mat');
%end
toc

K=[7 5 3];
K=[3 3 3];

N=[11 11 11 3];

y=y(1:N(1),1:N(2),1:N(3),1:N(4));
%z=hankel(y,K);
%visReconstruction(log(z(:,:)))

ND=3;
%HANKEL INDEXES
iw=1:prod(K);
sw=ind2subV(K,iw);
iX=cell(1,ND);
%parE.NC
for s=1:ND
    iX{s}=0:N(s)-K(s);
    iX{s}=bsxfun(@plus,sw(:,s),iX{s});
end
iXi=cell(1,ND);


NSa=prod(N(1:ND)-K(1:ND)+1);
An=zeros([NSa prod(K) N(4)],'like',y);
for n=1:prod(K)
    for s=1:ND;iXi{s}=iX{s}(n,:);end
    %if n==1
    %    iXi{1}(1)
    %    iXi{2}(1)
    %    iXi{3}(1)
    %end
    An=dynInd(An,n,2,reshape(dynInd(y,iXi,1:ND),[NSa 1 N(4)]));
end
visReconstruction(log(An(:,:)))
%return


NZZ=1;
if ND==3;NZZ=N(3)-K(3)+1;end
NSa=prod(N(1:2)-K(1:2)+1);
A=[];

for zz=1:NZZ
    Ab=zeros([NSa prod(K) N(4)],'like',y);
    for w=1:prod(K)
        for s=1:2;iXi{s}=iX{s}(w,:);end
        iXi{3}=iX{3}(w,zz);
        %if w==1 && zz==1
        %    iXi{1}(1)
        %    iXi{2}(1)
        %    iXi{3}(1)
        %    pause
        %end
        Ab=dynInd(Ab,w,2,reshape(dynInd(y,iXi,1:3),[NSa 1 N(4)]));
    end
    A=cat(1,A,Ab);
end
visReconstruction(log(A(:,:)))