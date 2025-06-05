classdef preluLayer < nnet.layer.Layer
    % Example custom PReLU layer.

    properties (Learnable)
        % Layer learnable parameters
            
        % Scaling coefficient
        Alpha
    end
    
    methods
        function layer = preluLayer(numChannels, name) 
            % layer = preluLayer(numChannels, name) creates a PReLU layer
            % with numChannels channels and specifies the layer name.

            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "PReLU with " + numChannels + " channels";
        
            % Initialize scaling coefficient.
            layer.Alpha = rand([1 1 numChannels]); 
        end
        
        function Z = predict(layer, X)
            % Z = predict(layer, X) forwards the input data X through the
            % layer and outputs the result Z.
            
            Z = max(0, X) + layer.Alpha .* min(0, X);
        end
        
        function [dLdX, dLdAlpha] = backward(layer, X, Z, dLdZ, memory)
            % [dLdX, dLdAlpha] = backward(layer, X, Z, dLdZ, memory)
            % backward propagates the derivative of the loss function
            % through the layer.
            % Inputs:
            %         layer    - Layer to backward propagate through
            %         X        - Input data
            %         Z        - Output of layer forward function
            %         dLdZ     - Gradient propagated from the deeper layer
            %         memory   - Memory value which can be used in backward
            %                    propagation
            % Outputs:
            %         dLdX     - Derivative of the loss with respect to the
            %                    input data
            %         dLdAlpha - Derivative of the loss with respect to the
            %                    learnable parameter Alpha
            
            dLdX = layer.Alpha .* dLdZ;
            dLdX(X>0) = dLdZ(X>0);
            dLdAlpha = min(0,X) .* dLdZ;
            dLdAlpha = sum(sum(dLdAlpha,1),2);
            
            % Sum over all observations in mini-batch.
            dLdAlpha = sum(dLdAlpha,4);
        end
    end
end
