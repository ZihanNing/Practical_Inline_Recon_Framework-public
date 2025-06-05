function [anom,maxProfRisk]=anomalyDetection(x,pathOu,path,name,typ,no,pda,sig,Assign,anom,dyn)

%ANOMALYDETECTION   Performs anomaly detection on spectra
%   [ANOM,MAXPROFRISK]=ANOMALYDETECTION(X,PATHOU,PATH,NAME,TYP,{NO},{PDA},{SIG},{ASSIGN},{ANOM},{DYN})
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
if nargin<11 || isempty(dyn);dyn=[0 0 0];end

path=strrep(path,'/','-');
gpu=isa(x,'gpuArray');
N=size(x);N(end+1:4)=1;

NXO=N;

%Tests to look for the shape of the spike on 2019_05_13/CO_99330, volume
%Par2Rec{11}=66
% xa=dynInd(x,[22 12],[2 8]);xa=squeeze(xa);
% xb=xa(1:41,:);xb=fftshift(fftGPU(ifftshift(xb,1),1),1);
% xa=log(abs(xa)).*sign(xa);
% xa=xa.';xb=xb.';
% 
% visReconstruction(xa,0,[],[],'Readout','Coils')
% visReconstruction(xa(:,1:41),0,[],[],'Zoomed readout','Coils')
% figure
% plot(abs(xa(23,1:41)))
% title('Log-magnitude')
% figure
% plot(exp(abs(xa(23,1:41))))
% title('Magnitude')
% figure
% plot(angle(xa(23,1:41)));
% title('Phase')
% figure
% plot(abs(xa(:,1:41))')
% title('Log-magnitude')
% figure
% plot(abs(xb(23,:))')
% title('Spatial magnitude')
% visReconstruction(xb,[],[],[],'Readout','Coils')
% 
% pause


if typ>=2
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
    
    if gpu;x=gpuArray(x);end
    
    %WE RESHAPE TO GET THE DIFFERENT NOISE PROFILES ACQUIRED
    N=size(x);N(end+1:4)=1;
    Nno=size(no);
    Nno(2)=Nno(1)/N(1);
    Nno(1)=N(1);
    no=reshape(no,Nno);

    %WE GET THE SPECTRAL NOISE AND SIGNAL RESPONSE
    xF=gather(x);
    no=sqrt(N(1))*ifftGPU(ifftshift(no,1),1);
    for n=1:N(4);x=dynInd(x,n,4,sqrt(N(1))*ifftGPU(ifftshift(dynInd(x,n,4),1),1));end    
    
    %WE NORMALIZE TO THE LEVEL IN THE FOV OF INTEREST
    rangeNoiseP=[0.45 0.95];%Top values discarding potential outliers    
    rangeNoise=floor(rangeNoiseP(1)*N(1))+1:floor(rangeNoiseP(2)*N(1));
    [~,indSort]=sort(multDimSum(abs(no).^2,[2 4]));
    noRobust=dynInd(no,indSort(rangeNoise),1);
    [~,noCov4]=standardizeCoils(noRobust,noRobust);

    %WE PRODUCE A STANDARDIZATION ALONG THE CHANNELS
    noStd4=standardizeCoils(no,noCov4);
 
    %WE CALCULATE THE RECEIVER RESPONSE    
    K=2*prod(Nno([2 4]));    
    noStd4Av24=multDimSum(abs(noStd4).^2,[2 4])/K;
    %figure
    %plot(noStd4Av24)
    %hold on
    
    
    p=0.005;%p-value
    p=p/N(1);%Bonferroni correction for multiple hypothesis
    th=2*gammaincinv(1-p,K/2)/K;%To force that they are recognized as outliers later on
    noStd4Av24(noStd4Av24>2*th)=2*th;%For extreme outliered noise profiles, such as in 2018_09_03/HE_43130/he_03092018_1422427_15_2_dhcp8angiocowsenseV4.raw    
    noStd4Av24F=cdfFilt(noStd4Av24,'med',17,'circular');
    noStd4Av24N=K*noStd4Av24./noStd4Av24F;
    noOutlProb=1-gammainc(noStd4Av24N/2,K/2);
    noOutl=noOutlProb<p;
    if sum(noOutl)>0;fprintf('OUTLIERS SEEN IN NOISE SCAN\n');end
    noStd4Av24O=noStd4Av24;
    noStd4Av24O(noOutl)=noStd4Av24F(noOutl);%Replacement of outliers
    H=buildFilter(N(1),'tukey',20/N(1),gpu,0.5);
    noCov1=real(filtering(noStd4Av24O,H));
    %figure
    %plot(noStd4Av24)
    %hold on
    %plot(noCov1)
    %pause
    %pause
    
       
    %WE COMPUTE THE NOISE ONLY FOV    
    xf=standardizeCoils(x,noCov4);
    xf=bsxfun(@rdivide,xf,sqrt(noCov1));
    xf4=sum(abs(xf).^2,4)/N(4);xf=[];    
    minRatBack=64;%64;%64;%64;%32;%Minimum ratio of assumed background area (as a percentage of the FOV)
    %outlRangePerc=[0.99 1];%[0.1 0.5];%Originally range where outliers are not expected but we've used the full range to escape from areas with signal
    outlRangePerc=0.001*N(2);%[0.1 0.5];%Originally range where outliers are not expected but we've used the full range to escape from areas with signal    
    outlRangePerc=max(ceil(outlRangePerc),8);  
    outlRange=N(2)-outlRangePerc+1:N(2);
    %fprintf('Range maximum spectra: %d\n',length(outlRange));
    xf4s=sort(xf4,2);%Sorted energies
    xf4se=xf4s(:,outlRange);%Profiles in the range
    xf4se2=sum(xf4se,2)/N(2);%Energy over the profiles
    NB=floor(N(1)/minRatBack);%Size of the background
    NB=max(NB,4);
    ratFilter=8;
    NF=floor(N(1)/ratFilter);%Size of window to estimate the background
    NF=max(NF,16);
    wD=ones([NF 1],'like',xf4);%Background window (to convolve)
    xf4se2=convnfft(xf4se2,wD);%Energy of rectangular areas with size N(1)/minRatBack
    [~,iminxc]=min(xf4se2,[],1);%Minimum energy over all possible windows of size NB
    wSl=[-ceil((NB-1)/2) floor((NB-1)/2)]+iminxc;%Window indexes
    wC=ones([3*N(1) 1],'like',xf4);%New window
    wCI=wC;%Indexes
    wCI(:)=(1:3*N(1))-N(1);%To prevent wraps
    wC(wCI<wSl(1) | wCI>wSl(2))=0;%Ouside the background 
    wC=reshape(wC,N(1),[]);%To place back into original space
    wC=sum(wC,2);%Sum over the replicas
    x1=dynInd(x,wC==1,1);
    %fprintf('Range background size: %d\n',sum(wC(:)));
    
    %WE GENERATE A SYNTHETIC NOISE SAMPLE
    pda=squeeze(pda);
    pdaun=unique(gather(pda),'rows');
    if gpu;pdaun=gpuArray(pdaun);end
    pda=permute(pda,[3 1 4 2]);
    pdaun=permute(pdaun,[1 3 4 2]);
    Nno=N;Nno(2)=5000;
    x=gather(x);
    
    nosF=randn(Nno,'like',real(x))+1i*randn(Nno,'like',real(x));
    nosF=bsxfun(@times,nosF,sqrt(noCov1));                
    nosF=standardizeCoils(nosF,noCov4,1);                
    nosF=fftGPU(nosF,1)/sqrt(N(1));
    for p=1:size(pdaun,1)    
        indPdaUn=bsxfun(@eq,pda,dynInd(pdaun,p,1));
        indPdaUn=sum(indPdaUn,4);
        
        nos=bsxfun(@rdivide,nosF,dynInd(pdaun,p,1));                          
        nos=round(nos);        
        nos=bsxfun(@times,nos,dynInd(pdaun,p,1));       
        nos=sqrt(N(1))*ifftGPU(nos,1);
        noRobust=dynInd(nos,wC==1,1);

        [~,noCov4P{p}]=standardizeCoils(noRobust,noRobust);
        nosStd4=standardizeCoils(nos,noCov4P{p});
        K=2*prod(Nno([2 4]));
        noCov1P{p}=multDimSum(abs(nosStd4).^2,[2 4])/K;
        noCov1P{p}=real(filtering(noCov1P{p},H));
        noCov1P{p}=dynInd(noCov1P{p},wC==1,1);

        %WE DECORRELATE INDEPENDENTLY FOR EACH PROFILE
        xp=dynInd(x1,indPdaUn~=0,2);
        xp=standardizeCoils(xp,noCov4P{p});
        xp=bsxfun(@rdivide,xp,sqrt(noCov1P{p}));
        x1=dynInd(x1,indPdaUn~=0,2,xp);
        nos=[];
    end
    nosF=[];    

    N=size(x1);N(end+1:4)=1;
    x12=abs(x1).^2;
    K=double(2*prod(N([1 4])));    
    x2s=double(multDimSum(x12,[1 4]));
    
    %STANDARDIZE-It does not hold for more than half the profiles outliered as in 2018_09_03/HE_43130/he_03092018_1410308_12_2_dhcp8mpragesenseV4.raw, so it is taken out
    x12=x12/(median(x2s)/K);
    
    %CHECKS
    %p=1e-12;    
    p=1e-9;
    x2s=double(multDimSum(x12,[1 4]));
    outlProbs=(1-gammainc(x2s/2,K/2))*N(2);
     
    %CORRECTION TO SUPPESS CENTER OF THE SPECTRUM IN VOLUMETRIC SCANS
    if size(Assign.z{3},14)>1
        kmax=2.5;
        k=sqrt(Assign.z{3}.^2+Assign.z{2}.^2);k=k(:)';
        NR=N(2)/length(k);
        k=repmat(k,[NR 1]);k=k(:)';
        outlProbs(k<kmax)=0.5*N(2);
    end
    maxProfRisk=gather(-log10(min(outlProbs)));    
    fprintf('Maximum profile risk: %0.2f\n',maxProfRisk);
    
    allSp=0;%To show the periodogram of all channels (1) or the spectrum and periodogram of all channels (2)
    
    if any(outlProbs<p)% || outlProb<p
        anom=anom+1;
        fprintf('OULIER DETECTED! Percentage of outlier profiles: %0.4f%%\n',100*sum(outlProbs<p)/N(2));       
        xW=[];
        if gpu;xF=gpuArray(xF);end%Spectrum
        if allSp~=2;xF=max(xF,[],4);end
        xF=log(abs(xF));
        xF=xF/max(xF(:));
        if allSp==2
            %xF=dynInd(xF,29,4);
            xF=permute(xF,[1 4 2 3]);
            xF=reshape(xF,[size(xF,1)*size(xF,2) size(xF,3),size(xF,4)]);
        end
        xF=repmat(xF,[1 1 3]);        
        xW=cat(1,xW,gather(xF));        
        xlim=xF(1:4,:,:);xlim(:)=0;
        xlim(:,:,3)=1;
        xW=cat(1,xW,gather(xlim));
        if gpu;x=gpuArray(x);end%Periodogram
        if allSp==0;x=max(x,[],4);end       
        x=log(abs(x));
        x=x/max(x(:));
        if allSp>0
            %x=dynInd(x,29,4);
            x=permute(x,[1 4 2 3]);
            x=reshape(x,[size(x,1)*size(x,2) size(x,3),size(x,4)]);
        end
        x=repmat(x,[1 1 3]);

        %xO=find(wC==1,1,'first')-1;if xO>=1;x(max(xO-1,1):xO,1:2:end,2)=1;end
        %xO=find(wC==1,1,'last')+1;if xO<=NXO(1);x(xO:min(xO+1:NXO(1)),1:2:end,2)=1;end

        %Seems wC is not working very well because it is sometimes
        %detecting the background at the center of the FOV, for instance
        %for diffusion tests on '2019_05_01/Co_97130';
        xO=find(abs(diff(wC))==1,1,'first');x(max(xO-1,1):xO,1:2:end,2)=1;
        xO=find(abs(diff(wC))==1,1,'last')+1;x(xO:min(xO+1:NXO(1)),1:2:end,2)=1;

        
        xW=cat(1,xW,gather(x)); 
        xW=cat(1,xW,gather(xlim));
        xF=double(outlProbs<p);
        xF=repmat(xF,[8*N(1) 1 3]);
        xF(:,:,2:3)=0;
        xW=cat(1,xW,gather(xF));
        xW=cat(1,xW,gather(xlim));      
        x1=sum(abs(x1).^2,1);  
        x1=repmat(x1,[2*N(1) 1 1 1]);
        x1=permute(x1,[1 4 2 3]);
        x1=reshape(x1,[],N(2));
        x1=log(x1);
        x1=x1/max(x1(:));
        x1=repmat(x1,[1 1 3]);        
        xW=cat(1,xW,gather(x1));  
        pathAnom=sprintf('%s/Anomalies',pathOu);
        if ~exist(pathAnom,'dir');mkdir(pathAnom);end
        imwrite(xW,sprintf('%s/%s-%s-%05d-%05d-%05d.png',pathAnom,path,name,dyn(1)+10000*allSp,dyn(2)+10000*allSp,dyn(3)+10000*allSp));   
        return;
    end
    %pause
elseif typ==1
    %WE EXTRACT POTENTIALLY ANOMALOUS POINTS
    xmax=max(abs(rec.z),[],4);
    xmax=bsxfun(@rdivide,xmax,median(xmax,1));
    xmax=bsxfun(@rdivide,xmax,median(xmax,2));
    xmaxmed=cdfFilt(xmax,'med',[3 1]);
    xmaxmed=max((xmax-xmaxmed)./(abs(xmaxmed)),0);

    % %1) WE CHECK FOR SPECIFICALLY DAMAGED CHANNELS
    % x=gather(x);
    % medx=median(x,4);
    % madx=median(abs(bsxfun(@minus,x,medx)),4);
    % outlSD=bsxfun(@rdivide,abs(bsxfun(@minus,x,medx)),madx+eps);%Stahel-Donoho outlyingness
    % outlSDN=outlSD;
    % outlSDN=bsxfun(@rdivide,outlSDN,multDimMed(outlSDN,1:3));%Normalized along the samples (assumption no more than half of the spectral samples are damaged for all damaged coils)
    % outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,4));%Normalized along the coils (assumption no more than half of the coils are damaged)
    % outlSDCh=outlSDN.^2;
    % outlSDCh=bsxfun(@times,outlSDCh,gather(xmaxmed));
    % outlSDCh=outlSDCh./prod(N(1:2));
    % outlSDCh=multDimSum(outlSDCh,1:2);
    % outlSDCh=outlSDCh./(median(outlSDCh)+eps);
    % thDamageCh=5;%It has broken at 2 althought I've introduced the normalization of xmax (lines 16 and 17 in the mean time)

    %2) WE CHECK FOR SPECIFICALLY DAMAGED READOUTS
    medx=median(xmax,1);
    madx=median(abs(bsxfun(@minus,xmax,medx)),1);
    outlSD=bsxfun(@rdivide,abs(bsxfun(@minus,xmax,medx)),madx+eps);%Stahel-Donoho outlyingness
    outlSDN=outlSD;
    outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,1));%Normalized along the phase encodes
    outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,2));%Normalized along the readouts (assumption no more than half of the readouts are damaged)
    outlSDRe=outlSDN.^2;
    outlSDRe=bsxfun(@times,outlSDRe,xmaxmed);
    outlSDRe=outlSDRe/prod(N(1:2));
    outlSDRe=multDimSum(outlSDRe,1);
    outlSDRe=outlSDRe./(median(outlSDRe)+eps);
    thDamageRe=2000;%It has broken at 1000

    %if any(outlSDRe>thDamageRe) || any(outlSDCh>thDamageCh)
    if any(outlSDRe>thDamageRe)
        anom=anom+1;
        fprintf('OULIER DETECTED\n');       
        x=max(x,[],4);
        x=log(abs(x));
        x=x/max(x(:));   
        pathAnom=sprintf('%s/Anomalies',pathOu);
        if ~exist(pathAnom,'dir');mkdir(pathAnom);end
        imwrite(gather(x),sprintf('%s/%s-%s-%d.png',pathAnom,path,name,anom));
        return; 
    end
end
