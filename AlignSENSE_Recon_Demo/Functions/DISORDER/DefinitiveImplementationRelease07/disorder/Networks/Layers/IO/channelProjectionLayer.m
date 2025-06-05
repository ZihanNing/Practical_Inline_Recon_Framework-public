classdef channelProjectionLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties.
        % Layer properties go here.
        
        %Dimensions of the array
        NChannels
    end

    %properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
    %end
    
    methods
        function layer = channelProjectionLayer(name,nChannels)
            % (Optional) Create a myLayer.
            % This function must have the same name as the layer.
            % Layer constructor function goes here.
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "Resampling of channels";
        
            % Set channels to be added size
            layer.NChannels = nChannels;           
            
        end
        
        function Z = predict(layer, X)
            % Forward input data through the layer at prediction time and
            % output the result.
            %
            % Inputs:
            %         layer    -    Layer to forward propagate through
            %         X        -    Input data
            % Output:
            %         Z        -    Output of layer forward function
            
            % Layer forward function for prediction goes here.
            
            NC=size(X);NC(end+1:4)=1;
            Z=reshape(X,[NC(1)/layer.NChannels layer.NChannels NC(2:3) NC(4)]);
            Z=permute(Z,[1 3 4 2 5]);
            Z=reshape(Z,[NC(1)/layer.NChannels NC(2) NC(3)*layer.NChannels NC(4)]);         
        end

        function [dLdX] = backward(layer, ~, ~, dLdZ, ~)
            % Backward propagate the derivative of the loss function through 
            % the layer.
            %
            % Inputs:
            %         layer             - Layer to backward propagate through
            %         X                 - Input data
            %         Z                 - Output of layer forward function            
            %         dLdZ              - Gradient propagated from the deeper layer
            %         memory            - Memory value from forward function
            % Outputs:
            %         dLdX              - Derivative of the loss with respect to the
            %                             input data
            %         dLdW1, ..., dLdWn - Derivatives of the loss with respect to each
            %                             learnable parameter
            
            % Layer backward function goes here.
                                    
            NC=size(dLdZ);NC(end+1:4)=1;
            dLdX=reshape(dLdZ,[NC(1) NC(2) NC(3)/layer.NChannels layer.NChannels NC(4)]);
            dLdX=permute(dLdX,[1 4 2 3 5]);
            dLdX=reshape(dLdX,[NC(1)*layer.NChannels NC(2) NC(3)/layer.NChannels NC(4)]);
        end
    end
end