classdef reshapingLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties.
        % Layer properties go here.
        
        %Dimensions of the array
        InputSize
        
        OutputSize
        
    end

    %properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
    %end
    
    methods
        function layer = reshapingLayer(name,inputSize,outputSize)
            % (Optional) Create a myLayer.
            % This function must have the same name as the layer.
            % Layer constructor function goes here.
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "Reshaping of dimensions";
        
            % Set array size
            layer.InputSize = inputSize;
            
            % Set array size
            layer.OutputSize = outputSize;
            
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
            
            NR=layer.OutputSize;NR(4)=size(X,4);        
            Z=reshape(X,NR);          
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
            
            NR=layer.InputSize;NR(4)=size(dLdZ,4);
            dLdX=reshape(dLdZ,NR);
        end
    end
end