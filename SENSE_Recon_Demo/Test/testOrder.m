addpath(genpath('/home/lcg13/Work/DISORDER7TR02'))

load('/home/lcg13/Work/DataDeepFetalR6/Par.mat','Par');

Lin=Par.Lin;Sli=Par.Sli;Seg=Par.Seg;slicePos=Par.slicePos;
minLin=min(Lin);minSli=min(Sli);
iminLin=find(Lin==minLin);iminSli=find(Sli==minSli);
SliList=Sli(iminLin);
slicePosList=slicePos(:,iminLin);
sliN=slicePosList(1:3,SliList==2)-slicePosList(1:3,SliList==1);
sliN=sliN/norm(sliN);
sliP=sum(slicePosList(1:3,:).*sliN,1);
[sliPmin,iSliPmin]=min(sliP);
iSliPmin=iminLin(iSliPmin);
sliP=sliP-sliPmin;
sliPmin2=mink(sliP,2);
fprintf('Slice separation alternative: %.2f\n',sliPmin2(2));
sliP=round(sliP/sliPmin2(2))+1;
Sli=sliP(Sli);%This is the actual order of slices
SliceOrder=Sli;
PEOrder=Lin;
NP=length(Sli);
NE=Par.NSeg;
fprintf('Number of echoes: %d\n',NE);
NK=length(unique(Lin));
NSl=length(unique(Sli));
fprintf('Number of k-space lines: %d\n',NK);
fprintf('Number of slices: %d\n',NSl);
PEsl=Lin(Sli==1);%Phase encodes for a given slice
SEsl=Seg(Sli==1);%Segments for a given slice
SHsl=zeros(1,NK);
p=[0 find(diff(SEsl)<0) NK];%Indexes of interfaces between shots
NSh=length(p)-1;
fprintf('Number of shots: %d\n',NSh);
for s=1:NSh;SHsl(p(s)+1:p(s+1))=s;end
ShFull=zeros(1,Par.NLin);
ShFull(PEsl)=SHsl;
AShots=zeros(Par.NLin,NSh);
for n=1:Par.NLin
    if ShFull(n)~=0;AShots(n,ShFull(n))=1;end
end
%figure
%imshow(AShots,[]);
for e=1:NE
    if sum(Seg==e)==NSh*NSl;break;end
end

ShotOrder=ShFull(Lin);
SliceOrderReduced=SliceOrder(Seg==e);
ShotOrderReduced=ShotOrder(Seg==e);

Time2ShotSlice=[ShotOrder' SliceOrder'];
Time2ShotSliceReduced=[ShotOrderReduced' SliceOrderReduced'];
ShotSlice2Time=zeros(NSh,NSl);
iShotSliceReduced=sub2indV([NSh NSl],Time2ShotSliceReduced);
ShotSlice2Time(iShotSliceReduced)=1:NSh*NSl;

Time2ShotSliceProcess=Time2ShotSlice;
NSS=NSh*NSl;
c=1;
while 1
    dt=Time2ShotSliceProcess(2,:)-Time2ShotSliceProcess(1,:);
    i=find(dt==0);
    transitionsF=diff([0;find(diff(Time2ShotSliceProcess(:,i)));NSS]);
    NSS=length(transitionsF);    
    if NSS==1;break;end
    groupSamples{c}=repelem(1:NSS,transitionsF');
    transitionsP=1+[0;find(diff(Time2ShotSliceProcess(:,i)))];
    Time2ShotSliceProcess=Time2ShotSliceProcess(transitionsP,:);
    c=c+1;
    if c==10;error('Maximum number of iterations reached without convergence of mapping scheme');end
end
    

figure
imshow(ShotSlice2Time,[])
ylabel('Shots','FontSize',16)
xlabel('Slices','FontSize',16)
colormap(jet)
c=colorbar;
c.Label.String='Time';
c.Label.FontSize=16;


%figure
%imshow(repmat(ShFull(ShFull~=0)',[1 16]),[])
%ylabel('Phase encoding axis','FontSize',16)
%colormap(jet)
%c=colorbar;
%c.Label.String='Shot';
%c.Label.FontSize=16;
