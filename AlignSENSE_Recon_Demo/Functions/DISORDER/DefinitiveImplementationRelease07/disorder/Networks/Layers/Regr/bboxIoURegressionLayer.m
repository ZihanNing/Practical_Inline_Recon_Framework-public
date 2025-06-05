classdef bboxIoURegressionLayer < nnet.layer.RegressionLayer
    % Example custom regression layer with mean-absolute-error loss.
    
    properties
        LogIoU

        % Layer properties go here.
    end
    
    methods
        function layer = bboxIoURegressionLayer(name,logIoU)
            % layer = ellipsoidOverlapRegressionLayer(name) creates a layer
            % of ellipsoid overlap for regression
			
            % Set layer name.
            layer.Name = name;
            
            %To use logIoU
            layer.LogIoU = logIoU;
            

            % Set layer description.
            layer.Description = 'Bbox intersection over union regression';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the overlap loss between
            % the predictions Y and the training targets T.
            
            lT=2*T(:,:,4:6,:);%Length of the bounding box, targets
            lY=2*Y(:,:,4:6,:);%Length of the bounding box, predictions
            
            vT=prod(lT,3);%Volume of the bounding box, targets
            vY=prod(lY,3);%Volume of the bounding box, predictions           
            
            xTl=T(:,:,1:3,:)-T(:,:,4:6,:);
            xTr=T(:,:,1:3,:)+T(:,:,4:6,:);
            xYl=Y(:,:,1:3,:)-Y(:,:,4:6,:);
            xYr=Y(:,:,1:3,:)+Y(:,:,4:6,:);
            
            lI=max(min(cat(2,xTr,xYr),[],2)-max(cat(2,xTl,xYl),[],2),0);%Volume of the intersection
            I=prod(lI,3);
            U=vT+vY-I;
            N=size(I,4);
            IoU=I./U;
            if layer.LogIoU
                IoU(IoU<1e-3)=1;%To prevent numerical issues                     
                loss=sum(-log(IoU),4)/N;%Mean over the minibatch     
            else
                IoU(IoU<1e-3)=0;%To prevent numerical issues          
                loss=sum(1-IoU,4)/N;
            end
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            
            lT=2*T(:,:,4:6,:);%Length of the bounding box, targets
            lY=2*Y(:,:,4:6,:);%Length of the bounding box, predictions
            
            vT=prod(lT,3);%Volume of the bounding box, targets
            vY=prod(lY,3);%Volume of the bounding box, predictions    
            
            xTl=T(:,:,1:3,:)-T(:,:,4:6,:);
            xTr=T(:,:,1:3,:)+T(:,:,4:6,:);
            xYl=Y(:,:,1:3,:)-Y(:,:,4:6,:);
            xYr=Y(:,:,1:3,:)+Y(:,:,4:6,:);
                        
            lI=max(min(cat(2,xTr,xYr),[],2)-max(cat(2,xTl,xYl),[],2),0);%Volume of the intersection
            I=prod(lI,3);
            U=vT+vY-I;
            IoU=I./U;
            N=size(I,4);
            
            %INTERSECTION DERIVATIVES
            IA=bsxfun(@rdivide,I,lI);%Areas
            IA=repmat(IA,[1 1 2 1]);
            multR=ones(1,1,6);
            multL=-ones(1,1,6);multL(4:6)=1;                   
            dIdY=bsxfun(@times,IA,bsxfun(@times,cat(3,xYl>xTl,xYl>=xTl),multL)+bsxfun(@times,cat(3,xYr<xTr,xYr<xTr),multR));
                       
            %UNION DERIVATIVES
            dvYdY=2*repmat(bsxfun(@rdivide,vY,lY),[1 1 2 1]);
            dvYdY(:,:,1:3,:)=0;%Only the raddi influence the union
            dUdY=dvYdY-dIdY;
            
            %LOSS DERIVATIVES
            if layer.LogIoU
                dLdY=(bsxfun(@rdivide,dUdY,U)-bsxfun(@rdivide,dIdY,I))/N;
            else
                dLdY=(bsxfun(@rdivide,bsxfun(@times,dUdY,I),U.^2)-bsxfun(@rdivide,dIdY,U))/N;
            end
            IoU=repmat(IoU,[1 1 6 1]);
            dLdY(IoU<1e-3)=0;
        end
    end
end
%Check UnitBox: An Advanced Object Detection Network