function anom=anomalyDetectionPCA(x,pathOu,path,name,typ)
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

%WE WORK WITH THE MEAN OVER THE CHANNELS
%x=min(x,[],4);
%x=mean(x,4);
%x=abs(x);
%x=bsxfun(@minus,x,mean(x,1));
x=max(x,[],4);

%WE IMPLEMENT THE METHODS IN FAST, PARAMETER FREEE OUTLIER IDENTIFICATION FOR ROBUST PCA
%NORM
%xn=bsxfun(@rdivide,x,sqrt(sum(x.*conj(x),1)));

%MINIMUM ACUTE ANGLE
acang=zeros([1 N(2)],'like',real(x));
for n=1:N(2)   
    xo=x(:,n);
    auxang=abs(bsxfun(@minus,x,xo));   
    auxang=auxang.^2;
    
    auxang=sqrt(sum(auxang,1))./norm(xo);    
    auxang(n)=[];
    acang(n)=min(auxang,[],2);
end

% 
% al=0.05;
% th=((4*sqrt(pi)*gamma((N(1)+1)/2)*log(1/(1-al/2)))/(N(2)^2*gamma(N(1)/2)))^(1/(N(1)-1));
% th
% 
% th=pi/2+sqrt(2)*erfcinv(2*(1-1/(2*N(2)^2*(N(2)-1))))/sqrt(N(1)-2);
% th

%th

figure
plot(sort(acang))
%hold on
%plot(1:N(2),th*ones(1,N(2)),'r')
pause

%if any(acang>th)
%    fprintf('OULIER DETECTED\n');        
%    x=log(abs(x));
%    x=x/max(x(:));               
%    %x=x(:,outlSDRe>thDamage);        
%    imwrite(gather(x),sprintf('%s/../../Snapshots/%s-%s.png',pathOu,path,name));
%    return  
%end
