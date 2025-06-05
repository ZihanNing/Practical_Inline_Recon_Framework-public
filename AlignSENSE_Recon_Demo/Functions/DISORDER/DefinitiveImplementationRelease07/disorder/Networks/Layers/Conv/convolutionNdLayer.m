 classdef convolutionNdLayer < nnet.layer.Layer
     
    properties
        % (Optional) Layer properties.

        % Layer properties go here.
        NW;
        
        NX;
        
        BoundaryCondition;
        
        InitType;
    end

    properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
        
         Bias;
         
         Weights;
        
    end
     
    methods
        function layer = convolutionNdLayer(name,nw,nx,bc,initType) 

            % Set layer name
            layer.Name = name;

            % Set layer description
            layer.Description = 'convolutionNdLayer';
            
            layer.NW=nw;%Size of the convolution kernel, as WidthxHeightx...xNInputChannelsxNOutputChannels
            
            layer.NX=nx;%Size of spatial dimensions of the data
            
            layer.InitType=initType;
            if layer.InitType==0
                layer.Weights=randn(layer.NW)*sqrt(2)/prod([nw(1:end-2) nw(end)]);
            elseif layer.InitType==1
                layer.Weights=(randn(layer.NW)+1)*sqrt(2)/prod([nw(1:end-2) nw(end)]);
            elseif layer.InitType==2
                layer.Weights=0.01*(randn(layer.NW)+1)*sqrt(2)/prod([nw(1:end-2) nw(end)]);
            end                
            
            layer.Bias=zeros([1 1 layer.NW(end)]);
            
            layer.BoundaryCondition=bc;%Boundary condition: one of 0,'circular','replicate','symmetric'

        end

        function Z = predict(layer,X)
            % Forward input data through the layer and output the result
            
            ND=length(layer.NX);
            NinCh=size(X,3);
            NouCh=layer.NW(ND+2);
            Nobs=size(X,4);
            
            X=reshape(X,[layer.NX NinCh Nobs]);
            X=im2convm(X,layer.NW(1:ND),1,layer.BoundaryCondition);
            X=reshape(X,[prod(layer.NX) NinCh Nobs prod(layer.NW(1:ND))]);
            X=permute(X,[2 4 1 3]);%Input channels / Size of convolution / Size of image / Number of images  
            X=reshape(X,NinCh*prod(layer.NW(1:ND)),[]);
            
            
            W=reshape(layer.Weights,[prod(layer.NW(1:ND)) layer.NW(ND+1) NouCh]);
            W=permute(W,[3 2 1]);
            W=reshape(W,NouCh,[]);
            
            Z=W*X;%Output channels---Size of imagexNumber of observations   
            
            Z=reshape(Z,[NouCh prod(layer.NX) Nobs]);
            Z=permute(Z,[2 1 3]);
            Z=reshape(Z,[layer.NX NouCh Nobs]);          
            B=reshape(layer.Bias,[ones(1,ND) NouCh]);
            Z=bsxfun(@plus,Z,B);
            
            Z=reshape(Z,[layer.NX(1) prod(layer.NX(2:end)) NouCh Nobs]);
        end
        
        function [Z,memory] = forward(layer,X)
            % Forward input data through the layer and output the result
            
            ND=length(layer.NX);
            NinCh=size(X,3);
            NouCh=layer.NW(ND+2);
            Nobs=size(X,4);
            
            X=reshape(X,[layer.NX NinCh Nobs]);
            X=im2convm(X,layer.NW(1:ND),1,layer.BoundaryCondition);
            X=reshape(X,[prod(layer.NX) NinCh Nobs prod(layer.NW(1:ND))]);
            X=permute(X,[2 4 1 3]);%Input channels / Size of convolution / Size of image / Number of images  
            X=reshape(X,NinCh*prod(layer.NW(1:ND)),[]);
            memory=X;
            
            
            W=reshape(layer.Weights,[prod(layer.NW(1:ND)) layer.NW(ND+1) NouCh]);
            W=permute(W,[3 2 1]);
            W=reshape(W,NouCh,[]);
            
            Z=W*X;%Output channels---Size of imagexNumber of observations   
            
            Z=reshape(Z,[NouCh prod(layer.NX) Nobs]);
            Z=permute(Z,[2 1 3]);
            Z=reshape(Z,[layer.NX NouCh Nobs]);          
            B=reshape(layer.Bias,[ones(1,ND) NouCh]);
            Z=bsxfun(@plus,Z,B);
            
            Z=reshape(Z,[layer.NX(1) prod(layer.NX(2:end)) NouCh Nobs]);
        end

        function [dLdX,dLdBias,dLdWeights] = backward(layer,X,~,dLdZ,memory)
            % Backward propagate the derivative of the loss function through 
            % the layer 
            
            ND=length(layer.NX);
            NinCh=size(X,3);
            NouCh=size(dLdZ,3);
            Nobs=size(dLdZ,4);
            
            %Derivative of the biases
            dLdBias=multDimSum(dLdZ,[1:2 4]);
            
            %Derivative of the weights
            dLdZr=permute(dLdZ,[3 1 2 4]);
            dLdZr=reshape(dLdZr,NouCh,[]);            
            dLdWeights=dLdZr*memory';%Output channels---Input channelsxSize of convolution
            dLdWeights=reshape(dLdWeights,[NouCh NinCh prod(layer.NW(1:ND))]);
            dLdWeights=permute(dLdWeights,[3 2 1]);
            dLdWeights=reshape(dLdWeights,layer.NW);
            
            %Input gradient
            W=reshape(layer.Weights,[prod(layer.NW(1:ND))*NinCh NouCh]);
            dLdX=W*dLdZr;%Convolution sizexNinput channels---Image sizexNobservations            
            dLdX=reshape(dLdX,[prod(layer.NW(1:ND)) NinCh prod(layer.NX) Nobs]);
            dLdX=permute(dLdX,[3 2 4 1]);
            dLdX=reshape(dLdX,[layer.NX NinCh Nobs prod(layer.NW(1:ND))]);
            dLdX=im2convm(dLdX,layer.NW(1:ND),0,layer.BoundaryCondition);
            dLdX=reshape(dLdX,[layer.NX(1) prod(layer.NX(2:end)) NinCh Nobs]);
        end
    end
end