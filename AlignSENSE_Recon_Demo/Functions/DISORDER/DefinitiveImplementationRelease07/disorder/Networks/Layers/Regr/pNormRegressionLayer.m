classdef pNormRegressionLayer < nnet.layer.RegressionLayer
    % Example custom regression layer with mean-absolute-error loss.
    
    properties
        Order

        % Layer properties go here.
    end
    
    methods
        function layer = pNormRegressionLayer(name,order)
            % layer = maeRegressionLayer(name) creates a
            % mean-absolute-error regression layer and specifies the layer
            % name.
			
            % Set layer name.
            layer.Name = name;
            
            %Set layer order
            layer.Order = order;

            % Set layer description.
            layer.Description = 'P-value norm for regression';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the MAE loss between
            % the predictions Y and the training targets T.

            % Calculate MAE.
            R = size(Y,3);
            %Y(:,:,1,:)=Y(:,:,1,:)*2;
            %T(:,:,1,:)=T(:,:,1,:)*2;
            meanAbsoluteError = sum(abs(Y-T).^layer.Order,3)/R;
    
            % Take mean over mini-batch.
            N = size(Y,4);
            loss = sum(meanAbsoluteError)/N;
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % Returns the derivatives of the MAE loss with respect to the predictions Y

            %Y(:,:,1,:)=Y(:,:,1,:)*2;
            %T(:,:,1,:)=T(:,:,1,:)*2;
            R = size(Y,3);
            N = size(Y,4);            
            dLdY = layer.Order*(Y-T).^(layer.Order-1)/(N*R);
        end
    end
end