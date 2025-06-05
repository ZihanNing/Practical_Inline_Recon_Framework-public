function xou=virtualBodyCoil(xin)

%xin=single(randn(256,128,256,32)+1i*randn(256,128,256,32));
%xin=gpuArray(xin);
%tic
gpu=isa(xin,'gpuArray');

NX=size(xin);NX(end+1:4)=1;

xou=zeros(NX(1:3),'like',xin);
wou=zeros(NX(1:3),'like',xin);

%sub=[32 32 32];
%NW=64;

%sub=[8 8 8];
%sub=[32 32 32];
%sub=[16 16 16];
NW=[64 64 64];
%NW=[64 64 64];
sub=NW/4;
%NW=[8 8 8];
NC=1;

if ~any(NW==1)
    for n=1:3  
        %w{n}=ones([NW(n) 1],'single');
        w{n}=single(gausswin(NW(n)));
        if gpu;w{n}=gpuArray(w{n});end
        ww{n}=circconvmtx(w{n},NX(n),1);
        wo{n}=circshift(circconvmtx(w{n},NX(n)),-(NW(n)-1));
        ww{n}=ww{n}(1:sub(n):end,:);
        wo{n}=wo{n}([1 1+sub(n)],:);    
        wo{n}=prod(wo{n},1);
        wo{n}=circshift(wo{n},-sub(n));
        wo{n}=wo{n}(1:NW(n));
    end
    wo{1}=permute(wo{1},[2 1]);wo{1}=single(wo{1}>0);
    wo{2}=wo{2};wo{2}=single(wo{2}>0);
    wo{3}=permute(wo{3},[1 3 2]);wo{3}=single(wo{3}>0);
    wr{1}=w{1};
    wr{2}=permute(w{2},[2 1]);
    wr{3}=permute(w{3},[2 3 1]);
    wwr=bsxfun(@times,bsxfun(@times,wr{1},wr{2}),wr{3});

    NP=[size(ww{1},1) size(ww{2},1) size(ww{3},1)];
    for n=1:NP(1)
        wf{1}=permute(ww{1}(n,:),[2 1]);    
        indnz{1}=find(wf{1}(:)~=0);
        if indnz{1}(1)==1 && indnz{1}(end)==NX(1);a=find(diff(indnz{1})>1);indnz{1}=circshift(indnz{1},-a);end
        for m=1:NP(2)      
            wf{2}=ww{2}(m,:);
            indnz{2}=find(wf{2}(:)~=0);
            if indnz{2}(1)==1 && indnz{2}(end)==NX(2);a=find(diff(indnz{2})>1);indnz{2}=circshift(indnz{2},-a);end
            for o=1:NP(3)
                wf{3}=permute(ww{3}(o,:),[1 3 2]);
                indnz{3}=find(wf{3}(:)~=0);
                if indnz{3}(1)==1 && indnz{3}(end)==NX(3);a=find(diff(indnz{3})>1);indnz{3}=circshift(indnz{3},-a);end

                xx=dynInd(xin,indnz,1:3);       
                xx=bsxfun(@times,wwr,xx);
                xx=compressCoils(xx,NC,[],[],[],1);              
                %NXX=size(xx);NXX(end+1:4)=1;
                %B=CNCGUnwrapping(xx);
                %xx=bsxfun(@times,xx,exp(1i*B));

                %THIS WOULD BE THE PHASE CORRECTION BUT IT DOES NOT SEEM TO WORK WELL
                %P=bsxfun(@times,xx,1./wwr);
                %P=xx;
                %if o>1;Pref{3}=circshift(Pref{3},[0 0 -sub(3)]);pu=3;
                %elseif m>1;Pref{2}=circshift(Pref{2},[0 -sub(2) 0]);pu=2;
                %elseif n>1;Pref{1}=circshift(Pref{1},[-sub(1) 0 0]);pu=1;
                %else pu=0;
                %end
                %if pu~=0
                %    PM=P;PM=bsxfun(@times,PM,sqrt(wo{pu}));
                %    PPM=Pref{pu};PPM=bsxfun(@times,PPM,sqrt(wo{pu}));
                %    xx=alignPhase(xx,PM,PPM);
                %end
                %P=bsxfun(@times,xx,1./wwr);     
                %%P=xx;
                %if m==1 && o==1;Pref{1}=P;end
                %if o==1;Pref{2}=P;end
                %Pref{3}=P;
            
                xou=dynInd(xou,indnz,1:3,dynInd(xou,indnz,1:3)+dynInd(xx,1,4));
                wou=dynInd(wou,indnz,1:3,dynInd(wou,indnz,1:3)+wwr);
            end
        end
    end
    xou=bsxfun(@times,xou,1./wou);
elseif sum(NW==1)==1
    d=find(NW==1);
    for n=1:NX(d)
        xx=dynInd(xin,n,d);
        xx=compressCoils(xx,NC,[],[],[],1);
        if n>1;xxref=alignPhase(xx,xx,xxref);else xxref=xx;end
        xou=dynInd(xou,n,d,dynInd(xxref,1,4));
    end
else
    d=find(NW==1);
    for n=1:NX(d(1))
        for m=1:NX(d(2))
            xx=dynInd(xin,[n m],d);
            xx=compressCoils(xx,NC,[],[],[],1);
            if m>1;xxref=alignPhase(xx,xx,xxref);else xxref=xx;end
            xou=dynInd(xou,[n m],d,dynInd(xxref,1,4));
        end
        xx=dynInd(xou,n,d(1));
        if n>1;xxrefb=alignPhase(xx,xx,xxrefb);else xxrefb=xx;end
        xou=dynInd(xou,n,d(1),xxrefb);        
    end    
end

end

function xx=alignPhase(xx,P,Pref)
    gpu=isa(xx,'gpuArray');
    NXX=size(xx);NXX(end+1:4)=1;

    P=reshape(P,[prod(NXX(1:3)) NXX(4)]);P=P.';                
    Pref=reshape(Pref,[prod(NXX(1:3)) NXX(4)]);Pref=Pref.';
    PR=P*Pref';
    [U,~,V]=svd(gather(PR));
    if gpu;U=gpuArray(U);V=gpuArray(V);end    
    PR=V*U';
    xx=reshape(xx,[prod(NXX(1:3)) NXX(4)]);          
    xx=(PR*xx.').';          
    xx=reshape(xx,NXX);
end

