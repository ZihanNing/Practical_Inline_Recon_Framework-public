classdef resamplingLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties.
        % Layer properties go here.
        
        %Dimensions of the array
        InputSize
        
        OutputSize
        
        MirrorDims
    end

    %properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
    %end
    
    methods
        function layer = resamplingLayer(name,inputSize,outputSize,mirrorDims)
            % (Optional) Create a myLayer.
            % This function must have the same name as the layer.
            % Layer constructor function goes here.
            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "Resampling of dimensions";
        
            % Set array size
            layer.InputSize = inputSize;layer.InputSize(end+1:2)=1;
            
            % Set array size
            layer.OutputSize = outputSize;layer.OutputSize(end+1:2)=1;
            
            % Set dimensions to be mirrored
            layer.MirrorDims = mirrorDims;layer.MirrorDims(end+1:2)=0;layer.MirrorDims(layer.MirrorDims~=0)=2;
            
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
            Nin=[layer.InputSize NC(3:4)];
            Z=reshape(X,Nin);
            Z=real(resampling(Z,layer.OutputSize,0,layer.MirrorDims));
            Z=reshape(Z,[layer.OutputSize(1) prod(layer.OutputSize(2:end)) NC(3:4)]);            
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
            Nou=[layer.OutputSize NC(3:4)];
            dLdX=reshape(dLdZ,Nou);
            dLdX=real(resampling(dLdX,layer.InputSize,0,layer.MirrorDims))/prod(layer.InputSize./layer.OutputSize);
            dLdX=reshape(dLdX,[layer.InputSize(1) prod(layer.InputSize(2:end)) NC(3:4)]);
        end
    end
end