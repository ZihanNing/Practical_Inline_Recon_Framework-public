function visEmpiricalSamplesOld(cov,esd,x,beta)

%VISEMPIRICALSAMPLES   Visualizes the empirical samples computed with the spectrode method versus the histogram obtained by simulations
%   VISEMPIRICALSAMPLES(COV,ESD,ENC,X,BETA)
%   * COV is the population covariance pdf
%   * ESD is the empirical sample distribution
%   * ENC is the encoding information
%   * X is the data used to generate the population covariance
%   * BETA is the random matrix shape
%

gpu=isa(x,'gpuArray');
ND=numDims(x);ND=max(ND,3);
NPE=size(cov.Enc.kRange{1},1);
if NPE<=2
    x=dynInd(x,ones(1,ND-1),[1 3:ND]);
else
    x=dynInd(x,ones(1,ND-2),3:ND);
end
N=size(x);
x=repmat(x,[1 1 NPE round(prod(N)/(NPE*beta))]);
x=plugNoise(x);
Enc=cov.Enc;

for p=1:NPE
    Enc.AcqSize=cov.Enc.AcqSize(p,:);
    for n=1:length(cov.Enc.kRange);Enc.kRange{n}=cov.Enc.kRange{n}(p,:);end
    x=dynInd(x,p,3,margosianFilter(dynInd(x,p,3),Enc));
end
N=size(x);
M=prod(N(3:4));
P=prod(N(1:2));
x=reshape(x,[prod(N(1:2)) prod(N(3:4))]);
x=x/sqrt(2*M);
[S,U,V]=svdm(gather(x));
U=[];V=[];
S=diag(S);
S=S.^2;
figure
histogram(S,ceil(P/10),'Normalization','pdf')
hold on
plot(esd.grid,esd.dens)
pause