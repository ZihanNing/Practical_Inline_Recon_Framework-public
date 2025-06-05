function anom=anomalyDetectionPrev(x,pathOu,path,name,typ)
% 
% %ANOMALYDETECTION   Performs anomaly detection on spectra
% %   INDSPIKE=ANOMALYDETECTION(X)
% %   * X is the array on which to detect spikes
% %   * ANOM is a flag to indicate whether an anomaly was detected
% %
% 

anom=0;
path=strrep(path,'/','-');
gpu=isa(x,'gpuArray');
N=size(x);N(end+1:4)=1;

%1) WE CHECK WHETHER THERE ARE SPECIFIC DAMAGED CHANNELS

%WE PERFORM A MEDIAN FILTER TO WORK ON CONFLICTING AREAS
xmax=max(abs(x),[],4);
% xmed=cdfFilt(xmax,'med',[3 1]);
% %xmed=abs(xmax-xmed)./abs(xmed);
% xmed=max((xmax-xmed)./abs(xmed),0);
% 
% x=gather(x);
% 
% %WE COMPUTE THE STAHEL-DONOHO OUTLYINGNESS
% medx=median(x,4);
% madx=median(abs(bsxfun(@minus,x,medx)),4);
% outlSD=bsxfun(@rdivide,abs(bsxfun(@minus,x,medx)),madx+eps);%Stahel-Donoho outlyingness
% 
% outlSDN=outlSD;
% outlSDN=bsxfun(@rdivide,outlSDN,multDimMed(outlSDN,1:3));%Normalized along the samples (assumption no more than half of the spectral samples are damaged for all damaged coils)
% outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,4));%Normalized along the coils (assumption no more than half of the coils are damaged)
% 
% outlSDCh=outlSDN.^2;
% outlSDCh=bsxfun(@times,outlSDCh,xmed);
% outlSDCh=outlSDCh./prod(N(1:2));
% outlSDCh=multDimSum(outlSDCh,1:2);
% outlSDCh=outlSDCh./(median(outlSDCh)+eps);
% 
% thDamage=5;%It has broken at 2
% badChannels=find(outlSDCh(:)>=thDamage);
% %goodChannels=find(outlSDCh(:)<thDamage);
% allCh=0;

%2) WE CHECK FOR SPIKES IN ALL CHANNELS
%xmed=abs(x);
%xmed=median(xmed,4);
xmed=bsxfun(@rdivide,xmed,median(xmed,1));
xmed=bsxfun(@rdivide,xmed,median(xmed,2));
xmedr=cdfFilt(xmed,'med',[3 1]);
xmedr=max((xmed-xmedr)./abs(xmedr),0);

%STAHEL-DONOHO OUTLYINGNESS OF EACH READOUT
medx=median(xmed,1);
madx=median(abs(bsxfun(@minus,xmed,medx)),1);
outlSD=bsxfun(@rdivide,abs(bsxfun(@minus,xmed,medx)),madx+eps);%Stahel-Donoho outlyingness

outlSDN=outlSD;
outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,1));%Normalized along the phase encodes
outlSDN=bsxfun(@rdivide,outlSDN,median(outlSDN,2));%Normalized along the readouts (assumption no more than half of the readouts are damaged)

outlSDRe=outlSDN.^2;
outlSDRe=bsxfun(@times,outlSDRe,xmedr);
outlSDRe=outlSDRe/prod(N(1:2));
outlSDRe=multDimSum(outlSDRe,1);
outlSDRe=outlSDRe./(median(outlSDRe)+eps);

%figure
%plot(log10(outlSDRe))
%pause

%figure
%imshow(log(abs(xmed)),[])
%pause

thDamage=2000;%It has broken at 2
if any(outlSDRe>thDamage)
    badChannels=1:length(outlSDCh);
    %goodChannels=[];
    allCh=1;
end

if ~isempty(badChannels)
    anom=1;
    if typ==2
        fprintf('OULIER AT CHANNELS:%s\n',sprintf(' %d',badChannels));
        xmax=gather(xmax);
        xmax=xmax/max(xmax(:));       
        imwrite(xmax,sprintf('%s/../../Snapshots/%s-%s-%d.png',pathOu,path,name,allCh));
        return;
    end
    
%     if ~isempty(goodChannels)
%         outlSDNgood=dynInd(outlSDN,goodChannels,4);
%         outlSDNgood=median(outlSDNgood,4);
% 
%         outlSDNbad=dynInd(outlSDN,badChannels,4);
%         outlSDNbad=median(outlSDNbad,4);
% 
%         outlSDNbad=outlSDNbad./(outlSDNgood+eps);
% 
%         thOutl=3;
%         outlSDNbad=outlSDNbad>thOutl;
% 
%         %WE SAVE INFORMATION TO FILE
%         NO=sum(outlSDNbad(:));
%         fprintf('\nBroken coils:%s\n',sprintf(' %d',badChannels));
%         fprintf('Number of outliers: %d\n',NO);
%         fprintf('Ratio of outliers: %.4f%%\n\n',100*sum(outlSDNbad(:))/numel(outlSDNbad));
%         %imwrite(outlSDNbad,sprintf('/home/lcg13/Data/AnomaliesNew/%s-%s-Outliers.png',path,name));
% 
%         %WE DETECT PROBLEMATIC CHANNELS
%         N=size(x);
%         xr=reshape(x,[prod(N(1:3)) N(4)]);   
%         xr=xr(outlSDNbad==1,:);
%         medx=median(xr,2);
%         madx=median(abs(bsxfun(@minus,xr,medx)),2);
%         outlSD=bsxfun(@rdivide,abs(bsxfun(@minus,xr,medx)),madx+eps);%Stahel-Donoho outlyingness
%         outlSD=mean(outlSD,1);
%         x=x(:,:,:,badChannels);
%         x=log(abs(x));
%         for n=1:length(badChannels)
%             xn=x(:,:,:,n);
%             xn=xn/max(xn(:));
%             imwrite(xn,sprintf('/home/lcg13/Data/AnomaliesNew/%s-%s-LogChannel%d-MeanOutlyingness%05.2f.png',path,name,badChannels(n),outlSD(badChannels(n))));
%         end
%     end
end

