function svr=svrSliceWeights(svr,non,metric)

%SVRSLICEWEIGHTS   Computes the reliability of each slice
%   SVR=SVRSLICEWEIGHTS(SVR)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   SVR is a svr structure containing different views (svr.y), spacings (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)

if nargin<3 || isempty(metric);metric='MS';end%Other options are FMS for filtered mean squares and NC for normalized correlation

svr=svrResiduals(svr);

de=0;%1e-3;%Regularization weight

V=cell(1,svr.NV);M=cell(1,svr.NV);Vi=cell(1,svr.NV);Mi=cell(1,svr.NV);
for v=1:svr.NV
    if ~isfield(svr,'EnViews') || svr.EnViews(v)==1    
        svr.E{v}=svr.E{v}-svr.y{v};
        
        %if non==-1;svr.E{v}=svr.E{v}-sum(svr.My{v}(:).*svr.E{v}(:))/sum(svr.My{v}(:));end%WE REMOVE THE DC COMPONENT DUE TO REGULARIZATION PERHAPS!
        if strcmp(metric,'FMS') && svr.ParSVR.fracOrd~=0;svr.E{v}=filtering(svr.E{v},svr.G{v},[],[],1,svr.FY{v},svr.FYH{v});end%%%%%NOTE WE'VE INTRODUCED THE QUICK FLAG TO ACCELERATE IT
        %svr.My{v}=max(svr.My{v},1e-3);
        
        %TAKE INTO ACCOUNT THAT THIS SHOULD BE E-Y, THE RESIDUALS, PROBABLY AVERAGED!!!                
        %V{v}=multDimMea(svr.My{v}.*abs(svr.E{v}).^2,svr.id(v,1:2));      
        %M{v}=multDimMea(svr.My{v},svr.id(v,1:2));
        V{v}=multDimMea(max(svr.My{v},1e-8).*abs(svr.E{v}).^2,svr.id(v,1:2));%CHANGED FROM 1E-3
        M{v}=multDimMea(max(svr.My{v},1e-8),svr.id(v,1:2));%CHANGED FROM 1E-3
        
        %V{v}=V{v}./M{v};%THIS WOULD BE USED IN CASE THE EXCITATION MASKS ARE DISABLED

        %MULTIPLICATION WITH THE EXCITATION MASKS
        Vi{v}=bsxfun(@times,V{v},svr.MPack{v});
        Mi{v}=bsxfun(@times,M{v},svr.MPack{v});
        V{v}(:)=0;
        M{v}(:)=0;
        for p=1:size(Vi{v},5)    
            Vp=dynInd(Vi{v},p,5);
            Mp=dynInd(Mi{v},p,5);        
            Vp=bsxfun(@times,Vp,svr.MSlices{v}{p});
            Mp=bsxfun(@times,Mp,svr.MSlices{v}{p});        
            Vp=mean(Vp,svr.id(v,3));
            Mp=mean(Mp,svr.id(v,3));
            Vp=Vp./Mp;
            Vp=bsxfun(@times,Vp,svr.MSlices{v}{p});
            Mp=bsxfun(@times,Mp,svr.MSlices{v}{p});
            V{v}=V{v}+sum(Vp,6);
            M{v}=M{v}+sum(Mp,6);            
        end       
        %V{v}=sqrt(max(de,V{v}));
        V{v}=sqrt(max(median(V{v}(:)),V{v}));
        %svr.W{v}=V{v}.^((non-2)/2);
        %V{v}=mean(V{v}.^((non-2)/2),svr.id(v,3));  
        svr.W{v}=V{v}.^(non-2);
        V{v}=mean(svr.W{v},svr.id(v,3));  
    end    
end

%for v=1:svr.NV    
%    svr.E{v}=max(svr.My{v},1e-8).*abs(svr.E{v}).^2;
%    svr.E{v}=cat(4,svr.E{v},max(svr.My{v},1e-8).*abs(svr.y{v}).^2);
%end
%svrWriteData(svr,sprintf('Step%02d-ResidualsSpatialCoordinates',svr.ParSVR.Step),svr.E,1,[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;           

%TO SEE THE RESIDUALS IN IMAGE SPACE
%y=svr.y;
%x=svr.x;
%if isfield(svr,'Ex');svr.x=svr.Ex;else svr.x(:)=0;end
%for v=1:svr.NV;svr.y{v}=abs(svr.E{v}).^2;end
%svr=svrCG(svr,10,[],[],7,2);
%svr.y=y;
%svr.Ex=svr.x;
%svr.x=x;
%svrWriteData(svr,sprintf('Step%02d-Residuals',svr.ParSVR.Step),svr.Ex,[],[],[],1);svr.ParSVR.Step=svr.ParSVR.Step+1;

if isfield(svr,'EnViews');V=V(svr.EnViews);end
V=mean(cat(4,V{:}),4);

for v=1:svr.NV
    if ~isfield(svr,'EnViews') || svr.EnViews(v)==1;svr.W{v}=bsxfun(@rdivide,svr.W{v},V);end
end

%figure
%plot(svr.W{1}(:))
%pause
