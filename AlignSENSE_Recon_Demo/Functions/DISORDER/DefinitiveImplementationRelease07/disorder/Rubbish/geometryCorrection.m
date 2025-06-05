function rec=geometryCorrection(rec,MR)

%MAPGEOMETRY   Maps the geometries of a set of datasets into the geometry of another dataset
%   Y=MAPGEOMETRY(Y,MT,REC)
%   * FIL is the filename of the dataset that needs to be mapped
%   * TYP are the data types to be mapped
%   * REC is a reconstruction structure. At this stage it may contain the naming information (rec.Names), the status of the reconstruction (.Fail), 
%   the .lab information (rec.Par), the fixed plan information (rec.Plan), the dynamic plan information (rec.Dyn), the data information
%   (rec.(rec.Plan.Types)), the information for correction (rec.Corr.(rec.Plan.Types)), the informaton for further sorting 
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   * REC is the modified reconstruction structure with the datasets mapped
%

gpu=rec.Dyn.GPU;

N=size(rec.x);N(end+1:3)=1;N=N(1:3);

%%%%NOTE THIS HAS TO BE GENERALIZED FOR SEVERAL PES
M=inv(rec.Par.Mine.AIso);

rGrid=generateGrid(N,gpu,N,zeros(1,3));
[dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rGrid{1}(:),rGrid{2}(:),rGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;
deGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');
if gpu;M=gpuArray(M);end
deGrid=M*deGrid;
%size(deGrid)
deeGrid=zeros([N 3],'like',rec.x);
%rec.Par.Mine.Isoc
for m=1:3;deeGrid(:,:,:,m)=reshape(deGrid(m,:),N);end

GeoCorrPars=rec.Par.Labels.GeoCorrPars;


gx_radius=GeoCorrPars(1);
gx_c_coeffs=zeros(1,192);
gx_s_coeffs=zeros(1,192);
gy_c_coeffs=zeros(1,192);
gy_s_coeffs=zeros(1,192);
gz_c_coeffs=zeros(1,192);
gz_s_coeffs=zeros(1,192);

gx_c_coeffs(1:168)=GeoCorrPars(170:337);
gy_s_coeffs(1:168)=GeoCorrPars(2:169);
aux=GeoCorrPars(338:end-1);
aux=repmat(aux,[12 1]);
aux(2:end,:)=0;
gz_c_coeffs=aux(:)';
%NOT SURE WHAT IT MEANS THE LAST GEOCORRPAR

C{1}=reshape(gx_c_coeffs,12,16)';
S{1}=reshape(gx_s_coeffs,12,16)';
C{2}=reshape(gy_c_coeffs,12,16)';
S{2}=reshape(gy_s_coeffs,12,16)';
C{3}=reshape(gz_c_coeffs,12,16)';
S{3}=reshape(gz_s_coeffs,12,16)';

tic
B=callSolidHarmonics(deeGrid,C,S,gx_radius);
toc
%min(deeGrid(:))
%max(deeGrid(:))
%min(B(:))
%max(B(:))

%rec.x=deeGrid;
%rec.x=cat(4,deeGrid,B);
%rec.x=B;

for n=1:3;dGrid{n}=B(:,:,:,n);end;
M=rec.Par.Mine.AIso;
deGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');
if gpu;M=gpuArray(M);end
deGrid=M*deGrid;
deeGrid=zeros([N 3],'like',rec.x);
%rec.Par.Mine.Isoc
for m=1:3;deeGrid(:,:,:,m)=reshape(deGrid(m,:),N);end





%tic
%tic
%rec.x=griddata(deeGrid(:,:,:,1),deeGrid(:,:,:,2),deeGrid(:,:,:,3),rec.x,B(:,:,:,1),B(:,:,:,2),B(:,:,:,3),'linear');
%toc
rec.x=interpn(rec.x,deeGrid(:,:,:,1),deeGrid(:,:,:,2),deeGrid(:,:,:,3),'linear',0);

%rec.x=interpn(deeGrid(:,:,:,1),deeGrid(:,:,:,2),deeGrid(:,:,:,3),rec.x,deeGrid(:,:,:,1),deeGrid(:,:,:,2),deeGrid(:,:,:,3),'linear',0);

%toc




%figure
%imshow(gx_c_coeffs~=0,[])
%figure
%imshow(gy_s_coeffs~=0,[])
%figure
%imshow(gz_coeffs~=0,[])
%pause



% gx_ref_radius=MR.Parameter.GetValue('RC_geom_corr_gx_ref_radius')
% gx_field_c_coeffs=MR.Parameter.GetValue('RC_geom_corr_gx_field_c_coeffs');
% gx_field_s_coeffs=MR.Parameter.GetValue('RC_geom_corr_gx_field_s_coeffs');
% gy_ref_radius=MR.Parameter.GetValue('RC_geom_corr_gy_ref_radius');
% gy_field_c_coeffs=MR.Parameter.GetValue('RC_geom_corr_gy_field_c_coeffs');
% gy_field_s_coeffs=MR.Parameter.GetValue('RC_geom_corr_gy_field_s_coeffs');
% gz_ref_radius=MR.Parameter.GetValue('RC_geom_corr_gz_ref_radius');
% gz_field_coeffs=MR.Parameter.GetValue('RC_geom_corr_gz_field_coeffs');







% 
% 
% 
% MTT=MTx\MTy;%EQUIVALENT TO MTT  
% Nsource=size(x);Nsource(end+1:3)=1;
% Ndestin=size(y);Ndestin(end+1:3)=1;
% 
% rdGrid=generateGrid(Ndestin,gpu,Ndestin,ones(1,3));
% [dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rdGrid{1}(:),rdGrid{2}(:),rdGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;
% destinGrid=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');dGrid{4}=[];
% if gpu;MTT=gpuArray(MTT);end
% 
% sdGrid=generateGrid(Nsource,gpu,Nsource,ones(1,3));
% [sGrid{1},sGrid{2},sGrid{3}]=ndgrid(sdGrid{1}(:),sdGrid{2}(:),sdGrid{3}(:));
% 
% destinGrid=MTT*destinGrid;
% for m=1:3
%     dGrid{m}=reshape(destinGrid(m,:),Ndestin);
%     dGrid{m}(dGrid{m}<min(sGrid{m}(:)))=min(sGrid{m}(:));
%     dGrid{m}(dGrid{m}>max(sGrid{m}(:)))=max(sGrid{m}(:));
% end
% x=interpn(sGrid{1},sGrid{2},sGrid{3},x,dGrid{1},dGrid{2},dGrid{3},'linear',0);
% 
% GeoCorrPars=rec.Par.Labels.GeoCorrPars;
% NG=length(GeoCorrPars);
% Spl=[2 (NG-18)/2 (NG-18)/2 16];
% Spl=cumsum(Spl);
% 
% 
% r0=GeoCorrPars(1:Spl(1));assert(sum(r0~=0)==1,'Not 1 (%d) non-zero radious parameter',sum(r0~=0));r0=r0(r0~=0);
% for n=1:2;Gl{n}=GeoCorrPars(Spl(n)+1:Spl(n+1));assert(sum(G{n}~=0)~=25,'Not 25 (%d) non-zero Gl{%d} parameters',sum(Gl{n}~=0),n);Gl{n}=Gl{n}(Gl{n}~=0);end
% Gt=GeoCorrPars(Spl(3)+1:Spl(4));assert(sum(Gt~=0)==8,'Not 8 (%d) non-zero Gt parameters',sum(Gt~=0));Gt=Gt(Gt~=0);
% %Orders of the Legendre Polinomials (as deducted from observed parameters...)
% M=9;N=15;
% L=size(rec.x);
% for n=1:2
%     G{n}=single(zeros(L(n)));
%     if gpu;G{n}=gpuArray(G{n});end
% end
% Gx=single(zeros(N));
% Gy=single(zeros(N));
% Gz=single(zeros(N));
% 
% for n=1:2:N
%     Gx
