function b=biasFieldEstimation(x,th)

%BIASFIELDESTIMATION   Estimates a bias field based on the method described
%in [1] MV Milchenko, OS Pianykh, JM Tyler, "The Fast Automatic Algorithm 
%for Correction of MR Bias Field," J Magn Res Imag, 24:891-900, 2006.
%   B=BIASFIELDESTIMATION(X,TH)
%   * X is the input volume
%   * TH is a gradient thresholding percentile
%   * B is the corrected volume
%

if nargin<2 || isempty(th);th=0.5;end

gpu=isa(x,'gpuArray');

N=size(x);N(end+1:3)=1;

H=buildFilter(N,'tukeyIso',ones(1,3),gpu,1);
x=filtering(x,H);
x=log(x);

[grix{2},grix{1},grix{3}]=gradient(x);

perm=[1 2 3;
      2 3 1;
      3 1 2];
thg=zeros(1,3);
grx=cell(1,3);
w=cell(1,3);
for n=1:3
    grix{n}=real(grix{n});
    grx{n}=4*grix{n};
    for m=1:2
        d=perm(n,m+1);
        grx{n}=grx{n}+dynInd(grix{n},[2:N(d) 1],d)+dynInd(grix{n},[N(d) 1:N(d)-1],d);
    end
    grx{n}=grx{n}/8;
    w{n}=sort(abs(grx{n}(:)));
    thg(n)=w{n}(ceil(prod(N)*th));
    w{n}=abs(grx{n})<thg(n);
end;grix=[];

%visReconstruction(grx{1});
%visReconstruction(grx{2});
%visReconstruction(grx{3});

powe=[];
order=3;
for n=0:order
    for m=0:order
        for o=0:order
            if o+m+n<=order && o+m+n>0;powe=[powe;o m n];end
        end
    end
end
NC=size(powe,1);
if gpu;powe=gpuArray(powe);end
powe=permute(powe,[3 2 4 1]);

xGrid=generateGrid(N,gpu,N);
xGridP=cell(1,3);

for n=1:3;xGridP{n}=bsxfun(@power,xGrid{n},dynInd(powe,n,2));end

xGridF=bsxfun(@times,bsxfun(@times,xGridP{1},xGridP{2}),xGridP{3});

xGridG=cell(1,3);
for n=1:3;xGridG{n}=bsxfun(@times,dynInd(powe,n,2),bsxfun(@rdivide,xGridF,xGrid{n}+eps));end

b=zeros([1 1 1 NC],'like',x);
A=zeros([1 1 1 NC NC],'like',x);
for n=1:3;
    A=A+multDimSum(bsxfun(@times,bsxfun(@times,w{n},xGridG{n}),permute(xGridG{n},[1 2 3 5 4])),1:3);
    b=b+multDimSum(bsxfun(@times,w{n}.*grx{n},xGridG{n}),1:3);
end
A=permute(A,[4 5 1 2 3]);
b=permute(b,[4 5 1 2 3]);

c=A\b;
c=permute(c,[2 3 4 1]);
b=real(exp(sum(bsxfun(@times,c,xGridF),4)));

