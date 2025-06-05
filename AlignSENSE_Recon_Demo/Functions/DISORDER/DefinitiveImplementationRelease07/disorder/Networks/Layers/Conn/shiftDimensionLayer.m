classdef shiftDimensionLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties.
        % Layer properties go here.
        
        %Number of dimensions of the array
        NumberOfDimensions
        
        %Dimensions of the array
        ArraySize
    end

    %properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
    %end
    
    methods
        function layer = shiftDimensionLayer(arraySize,name)
            % (Optional) Create a myLayer.
            % This function must have the same name as the layer.
            % Layer constructor function goes here.
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "Shifting of dimensions";
        
            % Set array size
            layer.ArraySize = arraySize;
            
            % Set number of dimensions
            layer.NumberOfDimensions=length(layer.ArraySize);
            
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
            N=[layer.ArraySize NC(3:4)];
            Z=reshape(X,N);
            perm=[2:layer.NumberOfDimensions 1 layer.NumberOfDimensions+1:layer.NumberOfDimensions+2];
            Z=permute(Z,perm);
            Z=reshape(Z,N(2),[],NC(3),NC(4));
            
        end

        %function [Z, memory] = forward(layer, X)
        %    % (Optional) Forward input data through the layer at training
        %    % time and output the result and a memory value.
        %    %
        %    % Inputs:
        %    %         layer  - Layer to forward propagate through
        %    %         X      - Input data
        %    % Outputs:
        %    %         Z      - Output of layer forward function
        %    %         memory - Memory value for backward propagation
        %
        %    % Layer forward function for training goes here.
        %end

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
            N=[layer.ArraySize NC(3:4)];
            perm=[2:layer.NumberOfDimensions 1 layer.NumberOfDimensions+1:layer.NumberOfDimensions+2];
            dLdX=reshape(dLdZ,N(perm));
            dLdX=ipermute(dLdX,perm);
            dLdX=reshape(dLdX,N(1),[],NC(3),NC(4));
        end
    end
end