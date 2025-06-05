function [anom,percReads,percSlices,percVolumes,coils]=detuningDetection(x,pathOu,path,name,typ,no,pda,sig,Assign,anom,percReads,percSlices,percVolumes,coils,dyn)

%DETUNINGDETECTION   Performs detection of detuned coils
%   [ANOM,MAXPROFRISK]=DETUNINGDETECTION(X,PATHOU,PATH,NAME,TYP,{NO},{PDA},{SIG},{ASSIGN},{ANOM},{DYN})
%   * X is the array on which to detect spikes
%   * PATHOU is the path to the reconstructions
%   * PATH is the study being reconstructed
%   * NAME is the name of the file being reconstructed 
%   * TYP is the type of anomaly detection, see variable 
%   Alg.DetectAnomalies in file reconAlgorithm.m
%   * {NO} are the noise samples
%   * {PDA} indicates gains of readouts
%   * {SIG} indicates polarities of readouts
%   * {ASSIGN} indicates the sampling order
%   * {ANOM} is the current status of the anomaly counter
%   * {DYN} indicates the explored dynamic
%   ** ANOM serves to indicate whether an anomaly was detected
%   ** MAXPROFRISK is the maximum profile risk detected
% 

if nargin<6;no=[];end
if nargin<7;pda=[];end
if nargin<8;sig=[];end
if nargin<9;Assign=[];end
if nargin<10 || isempty(anom);anom=0;end
if nargin<11 || isempty(percReads);percReads=zeros(2,1);end
if nargin<12 || isempty(percSlices);percSlices=zeros(2,1);end
if nargin<13 || isempty(percVolumes);percVolumes=zeros(2,1);end
if nargin<14;coils=[];end
if nargin<15 || isempty(dyn);dyn=[0 0 0];end

path=strrep(path,'/','-');
gpu=isa(x,'gpuArray');
N=size(x);N(end+1:4)=1;
Nor=size(x);Nor(end+1:14)=1;
Nor=[prod(Nor(2:3)) Nor(8) prod(Nor([5:6 7:14]))];%Phase encodes, Slices, Volumes

x=gather(x);
x=resPop(x,2:3,prod(N(2:3)),2);
pda=resPop(pda,[2 3],prod(N(2:3)),2);
sig=resPop(sig,[2 3],prod(N(2:3)),2);
N=size(x);N(end+1:14)=1;
x=resPop(x,5:14,prod(N(5:14)),3);
pda=resPop(pda,5:14,prod(N(5:14)),3);
sig=resPop(sig,5:14,prod(N(5:14)),3);
N=size(x);N(end+1:14)=1;
x=resPop(x,2:3,prod(N(2:3)),2);
pda=resPop(pda,2:3,prod(N(2:3)),2);
sig=resPop(sig,2:3,prod(N(2:3)),2);
sig=dynInd(sig,1,4);
x=dynInd(x,sig==-1,2,flip(dynInd(x,sig==-1,2),1));sig=[];        
xF=gather(x);

if gpu;x=gpuArray(x);end        
xr=sum(log(abs(x)+1),1);x=[];
%xr=bsxfun(@rdivide,xr,median(xr,4));
xr=bsxfun(@rdivide,xr,max(xr,[],4));

%xr=bsxfun(@rdivide,xr,median(xr,2));
xr=bsxfun(@rdivide,xr,max(xr,[],2));
%visReconstruction(squeeze(xr))

xr=1./xr;
%fprintf('Maximum value: %.4f\n',max(xr(:)));

p=2;
%[~,ir]=max(xr(:));
N=size(xr);
%ir=ind2subV(N,ir);ir=ir(4);
[xr,ir]=max(xr,[],4);

ir=unique(ir(xr>p));
coils=[coils ir(~ismember(ir,coils))];

allSp=3;
maxProfRisk=gather(max(xr(:)/p));
xrr=xr;
percReads=percReads+gather([sum(xrr>p);N(2)]);
xrr=max(reshape(xrr,Nor(1),[]),[],1);
percSlices=percSlices+gather([sum(xrr>p);N(2)/Nor(1)]);
xrr=max(reshape(xrr,Nor(2),[]),[],1);
percVolumes=percVolumes+gather([sum(xrr>p);N(2)/prod(Nor(1:2))]);
%fprintf('Percentage of detuned volumes: %.2f\n',100*percVolumes(1)/percVolumes(2));
%fprintf('Percentage of detuned slices: %.2f\n',100*percSlices(1)/percSlices(2));
%fprintf('Percentage of detuned readouts: %.2f\n',100*percReads(1)/percReads(2));

if any(xr>p)
    anom=anom+1;
    %fprintf('OULIER DETECTED! Percentage of outlier profiles: %0.4f%%\n',100*sum(xr>p)/N(2));       
    xW=[];
    if gpu;xF=gpuArray(xF);end%Spectrum
    if allSp==3;xF=dynInd(xF,ir,4);end
    xF=log(abs(xF));
    xF=xF/max(xF(:));
    xF=permute(xF,[1 4 2 3]);
    xF=reshape(xF,[size(xF,1)*size(xF,2) size(xF,3) size(xF,4)]);

    xF=repmat(xF,[1 1 3]);   
    xF(:,1:Nor(1):end,2)=1;
    xF(:,1:Nor(1):end,[1 3])=0;
    xF(:,Nor(1):Nor(1):end,2)=1;
    xF(:,Nor(1):Nor(1):end,[1 3])=0;
    
    xW=cat(1,xW,gather(xF));    
    xlim=xF(1:4,:,:);xlim(:)=0;
    xlim(:,:,3)=1;
    xW=cat(1,xW,gather(xlim));    
    xF=double(xr>p);
    xF=repmat(xF,[8*N(1) 1 3]);
    xF(:,:,2:3)=0;
    xW=cat(1,xW,gather(xF));
    xW=cat(1,xW,gather(xlim));  
    pathAnom=sprintf('%s/Detuning',pathOu);
    if ~exist(pathAnom,'dir');mkdir(pathAnom);end
    imwrite(xW,sprintf('%s/%s-%s-%05d-%05d-%05d%s.png',pathAnom,path,name,dyn(1)+10000*allSp,dyn(2)+10000*allSp,dyn(3)+10000*allSp,sprintf('-%02d',ir)));
end
