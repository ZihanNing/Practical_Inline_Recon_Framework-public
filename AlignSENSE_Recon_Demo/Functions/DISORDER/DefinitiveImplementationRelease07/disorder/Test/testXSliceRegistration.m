addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug07'));

load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/testXSliceRegistration');
ndT=6;
MB=3;
NT=size(T);

z=y;
vM=[NY(ND)-1:NY(ND) 1:7];
for m=vM
    tic
    z=dynInd(z,m,ND,solveXSliceRegistration(x,dynInd(y,m,ND),dynInd(T,m,ndT-2),Mk,Hsm,lambda,nX,toler));
    toc
end


Mk=eye(NY(3),'like',Mk);
permM=1:ndT;permM(1:2)=[3 ndT-1];permM([3 ndT-1])=1:2;
Mk=permute(Mk,permM);
w=y;
for m=vM
    tic
    for l=1:MB
        slEff=(1:NT(ndT-1))+(l-1)*NT(ndT-1);
        mEff=m+(l-1);if mEff>NY(ND);mEff=mEff-NY(ND);end
        if l==1;ym=dynInd(y,{slEff,mEff},[3 ND]);else ym=cat(3,ym,dynInd(y,{slEff,mEff},[3 ND]));end
        if l==1;Tm=dynInd(T,mEff,ndT-2);else Tm=cat(5,Tm,dynInd(T,mEff,ndT-2));end
    end
    ym=solveXSliceRegistration(x,ym,Tm,Mk,Hsm,lambda,nX,toler);
    for l=1:MB
        slEff=(1:NT(ndT-1))+(l-1)*NT(ndT-1);
        mEff=m+(l-1);if mEff>NY(ND);mEff=mEff-NY(ND);end
        w=dynInd(w,{slEff,mEff},[3 ND],dynInd(ym,slEff,3));
    end
    toc
end

r=y;
for m=vM
    tic
    for l=1:MB
        slEff=(1:NT(ndT-1))+(l-1)*NT(ndT-1);
        mEff=m+(l-1);if mEff>NY(ND);mEff=mEff-NY(ND);end
        if l==1;ym=dynInd(y,{slEff,mEff},[3 ND]);else ym=cat(3,ym,dynInd(y,{slEff,mEff},[3 ND]));end
        if l==1;Tm=dynInd(T,mEff,ndT-2);else Tm=cat(5,Tm,dynInd(T,mEff,ndT-2));end
    end
    ym=solveXSliceRegistration(x,ym,Tm,Mk,Hsm,lambda,nX,toler,21);
    for l=1:MB
        slEff=(1:NT(ndT-1))+(l-1)*NT(ndT-1);
        mEff=m+(l-1);if mEff>NY(ND);mEff=mEff-NY(ND);end
        r=dynInd(r,{slEff,mEff},[3 ND],dynInd(ym,slEff,3));
    end
    toc
end

for d=1:5
    close all
    visSegment(z,[],0,1,[],d);
    visSegment(w,[],0,1,[],d);
    visSegment(r,[],0,1,[],d);
    pause
end