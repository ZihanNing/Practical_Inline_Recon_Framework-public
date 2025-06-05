classdef diceBinaryPixelClassificationLayer < nnet.layer.ClassificationLayer
    % This layer implements the generalized dice loss function for training
    % semantic segmentation networks.
    
    properties(Constant)
        % Small constant to prevent division by zero. 
        Epsilon = 1e-8;
    end
    
    methods
        
        function layer = diceBinaryPixelClassificationLayer(name)
            % layer =  dicePixelClassificationLayer(name) creates a Dice
            % pixel classification layer with the specified name.
            
            % Set layer name.          
            layer.Name = name;
            
            % Set layer description.
            layer.Description = 'Dice loss';
        end
        
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the Dice loss between
            % the predictions Y and the training targets T.   

            % Weights by inverse of region size.
            W = 1 ./ (sum(sum(T,1),2).^2+layer.Epsilon);
            W=W-mean(W,3);
            W(W>0)=1;
            W(W==0)=0.5;
            W(W<0)=0;
            
            intersection = sum(sum(Y.*T,1),2);
            union = sum(sum(Y.^2 + T.^2, 1),2);          
            
            numer = 2*sum(W.*intersection,3) + layer.Epsilon;
            denom = sum(W.*union,3) + layer.Epsilon;
            
            % Compute Dice score.
            %dice = (numer./denom);
            dice = (numer./denom).^2;
            
            % Return average Dice loss.
            N = size(Y,4);
            loss = sum((1-dice))/N;
            
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % dLdY = backwardLoss(layer, Y, T) returns the derivatives of
            % the Dice loss with respect to the predictions Y.
            
            % Weights by inverse of region size.
            W = 1 ./ (sum(sum(T,1),2).^2+layer.Epsilon);
            W=W-mean(W,3);
            W(W>0)=1;
            W(W==0)=0.5;
            W(W<0)=0;
            
            intersection = sum(sum(Y.*T,1),2);
            union = sum(sum(Y.^2 + T.^2, 1),2);
     
            numer = 2*sum(W.*intersection,3) + layer.Epsilon;
            denom = sum(W.*union,3) + layer.Epsilon;
            
            N = size(Y,4);
      
            dLdY = (2*W.*Y.*numer./denom.^2 - 2*W.*T./denom)./N;
            
            dice=(numer./denom);
            dLdY=2*(1-dice).*dLdY;
        end
    end
end