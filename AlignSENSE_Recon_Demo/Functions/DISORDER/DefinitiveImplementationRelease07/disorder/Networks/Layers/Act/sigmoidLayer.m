 classdef sigmoidLayer < nnet.layer.Layer
    methods
        function layer = sigmoidLayer(name) 

            % Set layer name
            layer.Name = name;

            % Set layer description
            layer.Description = 'sigmoidLayer'; 

        end

        function Z = predict(~,X)
            % Forward input data through the layer and output the result

            Z = exp(X)./(exp(X)+1);
        end

        function dLdX = backward(~,~, Z,dLdZ,~)
            % Backward propagate the derivative of the loss function through 
            % the layer 

            dLdX = Z.*(1-Z) .* dLdZ;

        end
    end
end