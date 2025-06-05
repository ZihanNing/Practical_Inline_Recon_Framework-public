function T=computeFrequencyStabilization(x,subPyr,mirr,useAbs,rec)

% COMPUTEFREQUENCYSTABILIZATION compues frequency stabilization shifts
%   T=COMPUTEFREQUENCYSTABILIZATIONSHIFTS(X,SUBPYR,{MIRR}) 
%   * X is an array that needs to be stabilized, 4th dimension are dynamics 
%   and fifth dimension are components
%   * SUBPYR is the subsampling factor, serves to accelerate
%   * {MIRR} serves to mirror in the dynamic dimension, defaults to all 0
%   * {USEABS} compares only the magnitude
%   ** T are the stabilization shifts
%

ND=numDims(x);ND=max(ND,4);
if nargin<3 || isempty(mirr);mirr=zeros(1,ND);end
if nargin<4 || isempty(useAbs);useAbs=0;end
if nargin<5;rec=[];end

gpu=isa(x,'gpuArray');

N=size(x);N(end+1:ND)=1;
NEqLines=N(2);
NL=length(subPyr);

%N(4)==2 we stabilize against a reference, otherwise we do groupwisee
if N(4)~=2;T{2}=zeros([ones(1,3) N(4)],'like',x);else T{2}=0;end
for l=1:NL
    %GLOBAL SEARCH
    if N(4)~=2        
        sR=subPyr(l)*ones(1,4);sR(5:ND)=1;
        sR(4)=32;%8
    else
        sR=subPyr(l)*ones(1,3);sR(4:ND)=1;
    end        

    NR=ceil(N./sR);
    sRR=N./NR;
    %if numel(x)>4e8;x=gather(x);end    
    xl=resampling(x,NR,[],mirr,0);
    if gpu;xl=gpuArray(xl);x=gpuArray(x);end
    if useAbs~=2;xl=abs(xl);end
    %if N(5)>1 && N(4)==2%THIS WOULD BE POSSIBLE IF THE CHANNELS WERE COMPRESSED, BUT THEY ARE NOT
    %    NN=ceil(N(5)/subPyr(l));
    %    xl=dynInd(xl,1:NN,5);
    %end
    xl=fftGPU(xl,2);  
    kGrid=generateGrid([1 NR(2)],gpu,2*pi,ceil(([1 NR(2)]+1)/2));
    kGrid=ifftshift(kGrid{2},2);
    sub=2*(4.^(0:4));
    ran=1./(4.^(1:5));
    if N(4)~=2;TR{2}=real(resampling(T{2},[1 1 1 NR(4)],[],mirr(1:4),0));else TR{2}=T{2};end
    for q=1:length(ran)
        if N(4)~=2;xM=mean(shifting(xl,TR),4);else xM=dynInd(xl,1,4);end
        rGrid=-ran(q)*NR(2):1/(sub(q)):NR(2)*ran(q);
        rGrid=permute(rGrid,[1 3 4 5 6 2]);
        rGrid=bsxfun(@plus,rGrid,TR{2});
        if N(4)~=2;xll=xl;else xll=dynInd(xl,2,4);end
        if useAbs==1 || useAbs==0
            HH=exp(-1i*bsxfun(@times,rGrid,kGrid));   
            xll=bsxfun(@times,xll,HH);
            if useAbs==1
                xll=abs(ifftGPU(xll,2));          
                xM=abs(ifftGPU(xM,2));
            end
        else
            TT{2}=rGrid;            
            xll=ifftGPU(xll,2);
            xll=applyFrequencyStabilization(xll,TT,rec,[],NEqLines);
            xM=ifftGPU(xM,2);
        end
        xll=normm(xll,xM,[1:3 5]);                   
        [~,shi]=min(xll,[],6);
        TR{2}=indDim(rGrid,shi,6);        
    end   
    %figure;plot(rGrid(:),xll(:))
    %pause
    if N(4)~=2;T{2}=real(resampling(TR{2}*sRR(2),[1 1 1 N(4)],[],mirr(1:4),0));else T{2}=TR{2}*sRR(2);end
end
if N(4)~=2;T{2}=bsxfun(@minus,T{2},T{2}(1));end
