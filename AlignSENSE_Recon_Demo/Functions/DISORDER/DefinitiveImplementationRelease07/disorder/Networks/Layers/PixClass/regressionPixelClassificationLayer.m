classdef regressionPixelClassificationLayer < nnet.layer.ClassificationLayer
    % Example custom regression layer with mean-absolute-error loss.
    
    methods
        function layer = regressionPixelClassificationLayer(name)
            % layer = maeRegressionLayer(name) creates a
            % mean-absolute-error regression layer and specifies the layer
            % name.
			
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = 'Regression pixel classification';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the MAE loss between
            % the predictions Y and the training targets T.

            % Calculate error
            R = size(Y,3);
            meanSquaredError = sum((Y-T).^2,3)/R;
            meanSquaredError=sum(sum(meanSquaredError,1),2);
    
            % Take mean over mini-batch.
            N = size(Y,4);
            loss = sum(meanSquaredError)/N;
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % Returns the derivatives of the MAE loss with respect to the predictions Y

            R = size(Y,3);
            N = size(Y,4);          
            dLdY = 2*(Y-T)/(N*R);
        end
    end
end