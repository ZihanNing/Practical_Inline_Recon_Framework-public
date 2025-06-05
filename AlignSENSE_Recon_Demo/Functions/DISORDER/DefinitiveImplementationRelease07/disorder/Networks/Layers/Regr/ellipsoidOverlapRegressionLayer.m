classdef ellipsoidOverlapRegressionLayer < nnet.layer.RegressionLayer
    % Example custom regression layer with mean-absolute-error loss.
    
    methods
        function layer = ellipsoidOverlapRegressionLayer(name)
            % layer = ellipsoidOverlapRegressionLayer(name) creates a layer
            % of ellipsoid overlap for regression
			
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = 'Ellipsoid overlap regression';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the overlap loss between
            % the predictions Y and the training targets T.

            %First we map the predictions to the sphere given by the
            %training targets
            N=size(Y,4);
            Y(:,:,1:3,:)=(Y(:,:,1:3,:)-T(:,:,1:3,:))./T(:,:,4:6,:);
            Y(:,:,4:6,:)=Y(:,:,4:6,:)./T(:,:,4:6,:)-1;
            loss=sum(abs(Y(:)).^2)/N;%Mean over the Minibatch            
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % Returns the derivatives of the loss with respect to the predictions Y
            N=size(Y,4);
            Y(:,:,1:3,:)=(Y(:,:,1:3,:)-T(:,:,1:3,:))./T(:,:,4:6,:);
            Y(:,:,4:6,:)=Y(:,:,4:6,:)./T(:,:,4:6,:)-1;

            dLdY=2*Y/N;
        end
    end
end
%Check UnitBox: An Advanced Object Detection Network