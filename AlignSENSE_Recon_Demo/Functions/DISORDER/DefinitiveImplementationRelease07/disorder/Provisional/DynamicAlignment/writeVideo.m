function writeVideo(rec)

% WRITEVIDEO writes a video with the fMRI preprocessing information
%   WRITEVIDEO
%   * REC is a reconstruction structure. At this stage it may contain the 
%   naming information (rec.Names), the status of the reconstruction 
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), and the data 
%   information (rec.(rec.Plan.Types))
%

dyn1='x';
if isfield(rec,'e');dyn2='e';else dyn2='v';end
dyn3='M';
M=dynInd(rec.(dyn3),13,4);

%WE EXTRACT THE ROI OF THE RECONSTRUCTED DATA
x=extractROI(abs(rec.(dyn1)),rec.Enc.ROI,1,1:3);

%GENERATE THE FILE NAME
pathModa=fullfile(rec.Names.pathOu,numbe2Modal(rec.Par.Mine.Modal));
if ~exist(pathModa,'dir');mkdir(pathModa);end
fileName=fullfile(pathModa,rec.Names.Name);
suff=strcat(rec.Plan.Suff,rec.Plan.SuffOu);

N=size(rec.(dyn2));N(end+1:4)=1;
%WE CUT THROUGH DIFFERENT PLANES
perm=1:3;
A=[];
for n=1:3    
    permB=[perm 4];
    [xo,eo,mo]=parUnaFun({x,abs(rec.(dyn2)),M},@permute,permB);
    [xo,eo,mo]=parUnaFun({xo,eo,mo},@gather);
    cen=rec.Par.Mine.sphcen0(perm);
    perm=circshift(perm,[0 1]);
    [xo,eo,mo]=parUnaFun({xo,eo,mo},@dynInd,cen(3),3);
    NN=size(xo);NN=NN(1:2);
    NF=max(N(1:3))*ones(1,2);    
    [xo,eo,mo]=parUnaFun({xo,eo,mo},@padarray,NF-NN,0,'post');                
    B=bwboundaries(mo,'noholes');
    indB=sub2indV(NF,B{1});
    eo=reshape(eo,[prod(NF) N(4)]);
    eo(indB,:)=-1;
    eo=reshape(eo,[NF 1 N(4)]);
    B=cat(1,xo,eo);
    A=cat(2,A,B); 
end

A=A/max(A(:));
for l=1:N(4)
    B=A(:,:,:,l);
    indB=find(B<0);
    NB=size(B);
    B=reshape(B,[prod(NB(1:2)) 1]);
    B=repmat(B,[1 3]);
    B(indB,3)=1;
    B(indB,1:2)=0;
    B=reshape(B,[NB(1:2) 3]);
    [SIf,cm]=rgb2ind(B,256);
    if l==1
        imwrite(SIf,cm,sprintf('%s_Vi%s.gif',fileName,suff),'gif', 'Delay',0.2,'Loop',inf);
    else
        imwrite(SIf,cm,sprintf('%s_Vi%s.gif',fileName,suff),'gif', 'WriteMode','append','Delay',0.2);
    end
end