function svr=svrSliceWeightsAdaptive(svr,non,metric)

%SVRSLICEWEIGHTS   Computes the reliability of each slice
%   SVR=SVRSLICEWEIGHTS(SVR)
%   SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%   SVR is a svr structure containing different views (svr.y), spacings 
%   (svr.MS), orientations (svr.MT) and reconstruction parameters (svr.rec)
%

if nargin<3 || isempty(metric);metric='MS';end%Other options are FMS for filtered mean squares and NC for normalized correlation

svr=svrResiduals(svr);
gpu=isa(svr.y{1},'gpuArray');

de=1e-3;%Regularization weight

slR=[];slW=[];
V=cell(1,svr.NV);M=cell(1,svr.NV);Vi=cell(1,svr.NV);Mi=cell(1,svr.NV);
for v=1:svr.NV
    if ~isfield(svr,'EnViews') || svr.EnViews(v)==1    
        svr.E{v}=svr.E{v}-svr.y{v};                
        %if non==-1;svr.E{v}=svr.E{v}-sum(svr.My{v}(:).*svr.E{v}(:))/sum(svr.My{v}(:));end%WE REMOVE THE DC COMPONENT DUE TO REGULARIZATION PERHAPS!
        if strcmp(metric,'FMS') && svr.ParSVR.fracOrd~=0;svr.E{v}=filtering(svr.E{v},svr.G{v},[],[],1,svr.FY{v},svr.FYH{v});end%%%%%NOTE WE'VE INTRODUCED THE QUICK FLAG TO ACCELERATE IT
        
        
        V{v}=multDimMea(svr.My{v}.*abs(svr.E{v}).^2,svr.id(v,1:2));      
        M{v}=multDimMea(svr.My{v},svr.id(v,1:2));
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
         %if non<-1         
         %    sZ=numel(V{v})/svr.MBFactor(v);
         %    r=V{v}(1:sZ);
         %    w=M{v}(1:sZ);
         %    slR=cat(2,slR,r(:)');
         %    slW=cat(2,slW,w(:)');           
         %end
    end
end
  %if non<-1
  %    medslR=median(slR); 
  %end

%  if non<0
%      medslR=median(slR);
%      slR=slR-medslR;
%      save inpRobLoss slR slW
%      if non==-3
%          [d,c]=shapeParameterEstimation(slR,slW);
%      elseif non==-4
%          [c,d]=robustLossEstimation(double(gather(slR)),double(gather(slW)));       
%      end
%      fprintf('Shape parameter: %.4f\n',d);
%      fprintf('Scale parameter: %.4f\n',c);
%  end
for v=1:svr.NV
    if ~isfield(svr,'EnViews') || svr.EnViews(v)==1            
%         if non<-1
%             svr.W{v}=weightsRobustLoss(abs(V{v}),non,medslR)*medslR^2;
             %figure
             %plot(svr.W{v}(:));
             %pause
             
%         elseif non==-1          
%             V{v}=max(c,abs(V{v}-medslR));
%             svr.W{v}=V{v}.^((d-2)/2);
%             V{v}=mean(V{v}.^((d-2)/2),svr.id(v,3));  
%         else
            V{v}=sqrt(max(de,V{v}));
            %svr.W{v}=V{v}.^((non-2)/2);
            %V{v}=mean(V{v}.^((non-2)/2),svr.id(v,3));
            V{v}=V{v}.^(non-2);
            svr.W{v}=V{v};
            V{v}=mean(V{v},svr.id(v,3));
 %        end
    end
end

%if non>=-1
    if isfield(svr,'EnViews');V=V(svr.EnViews);end
    Vme=mean(cat(4,V{:}),4);
    for v=1:svr.NV        
        if ~isfield(svr,'EnViews') || svr.EnViews(v)==1        
            svr.W{v}=bsxfun(@rdivide,svr.W{v},Vme);
             %figure
             %plot(svr.W{v}(:));
             %pause
        end
    end
%end

function y=weightsRobustLoss(x,a,c,e)
    if nargin<4 || isempty(e);e=1e-5;end    
    
    b=abs(2-a)+e;
    d=a+e-2*e*(a<0);
    y=(c^(-2))*(((x/c).^2/b+1)).^(d/2-1);
end

end

