classdef waveletObj
% Class for output of a wavelet transform performed by the dwt class.
%
% Inspired by WAVELET class written by Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne, 17-07-2009
%
% (c) Corey Baron 2016

	properties
        high = [];   % high pass filter coeffs. NlxNd cell array, where Nl = # levels, Nd = all permutations over dimensions
        low = [];    % low pass filter coeffs
        extra = [];  % optional field for another transform etc
    end

	methods
		function obj = waveletObj(wavIn,wavTemplate) 
            if ~isstruct(wavIn) && isscalar(wavIn)
                % For a scalar input, just replicate the scalar into all
                % fields based on template
                obj = wavTemplate;
                obj.low = wavIn;
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = wavIn;
                    end
                end
                if isfield(obj,'extra')
                    obj.extra = wavIn;
                end
            else
                obj.high = wavIn.high;
                obj.low = wavIn.low;
                if isfield(wavIn,'extra')
                    obj.extra = wavIn.extra;
                end
            end
        end
        
        function obj = abs(obj)
            obj.low = abs(obj.low);
            for i = 1:length(obj.high)
                for m = 1:length(obj.high{i})
                    obj.high{i}{m} = abs(obj.high{i}{m});
                end
            end
            if ~isempty(obj.extra)
                obj.extra = abs(obj.extra);
            end
        end
        
        function obj = sign(obj)
            obj.low = sign(obj.low);
            for i = 1:length(obj.high)
                for m = 1:length(obj.high{i})
                    obj.high{i}{m} = sign(obj.high{i}{m});
                end
            end
            if ~isempty(obj.extra)
                obj.extra = sign(obj.extra);
            end
        end
        
        function obj2 = gather(obj)
            obj2.low = gather(obj.low);
            for i = 1:length(obj.high)
                for m = 1:length(obj.high{i})
                    obj2.high{i}{m} = gather(obj.high{i}{m});
                end
            end
            if ~isempty(obj.extra)
                obj2.extra = gather(obj.extra);
            else
                obj2.extra = [];
            end
        end
        
        function obj = uminus(obj) % -w
            obj.low = -obj.low;
            for i = 1:length(obj.high)
                for m = 1:length(obj.high{i})
                    obj.high{i}{m} = -obj.high{i}{m};
                end
            end
            obj.extra = -obj.extra;
        end

        function obj = mtimes(alpha,obj) 
            % Note: for different alphas at different levels/subbands,
            % alpha should be of the waveletObj class. Can have a different
            % scalar at each level or subband (and in extra field)
            if ~isa(alpha, 'waveletObj') && isa(obj, 'waveletObj')
                % case scalar * wavelet
                obj.low = alpha * obj.low;
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = alpha * obj.high{i}{m};
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = alpha * obj.extra;
                end
            elseif isa(alpha, 'waveletObj') && ~isa(obj, 'waveletObj') 
                % case wavelet * scalar
                obj = mtimes(obj,alpha);
            elseif isa(alpha, 'waveletObj') && isa(obj, 'waveletObj') 
                % case wavelet * wavelet
                obj.low = alpha.low .* obj.low;
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = alpha.high{i}{m} .* obj.high{i}{m};
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = alpha.extra .* obj.extra;
                end
            end
        end
        
        function obj = times(alpha,obj) 
            obj = mtimes(alpha,obj);
        end
        
        function obj = minus(alpha,obj)
            obj = -obj;
            obj = alpha + obj;
        end
        
        function obj = plus(alpha,obj)
            if ~isa(alpha, 'waveletObj') && isa(obj, 'waveletObj')
                % case scalar * wavelet
                obj.low = plus(alpha, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = plus(alpha, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = plus(alpha, obj.extra);
                end
            elseif isa(alpha, 'waveletObj') && ~isa(obj, 'waveletObj') 
                % case wavelet * scalar
                obj = plus(obj,alpha);
            elseif isa(alpha, 'waveletObj') && isa(obj, 'waveletObj') 
                % case wavelet * wavelet
                obj.low = plus(alpha.low, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = plus(alpha.high{i}{m}, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = plus(alpha.extra, obj.extra);
                end
            end
        end

        function out = subsref(obj,S)
            if isequal(S.type, '()') && isequal(S.subs,{':'}) % obj.(:)
                % Extract all the values
                out = obj.low(:);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        out = cat(1, out, obj.high{i}{m}(:));
                    end
                end
                if ~isempty(obj.extra)
                    out = cat(1, out, obj.extra(:));
                end
            elseif (length(S)==2) && isequal(S(1).type , '.') && isequal(S(1).subs , 'high')...
                    && isequal(S(2).type , '()') && isequal(S(2).subs,{':'}) % obj.high(:)
                % Extract all the high frequency values
                out = [];
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        out = cat(1, out, obj.high{i}{m}(:));
                    end
                end
            elseif (length(S)==3) && isequal(S(1).type , '.') && isequal(S(1).subs , 'high')...
                    && isequal(S(2).type , '{}') && (length(S(2).subs)==1)  && (length(S(2).subs{1})==1)...
                    && isequal(S(3).type , '()') && isequal(S(3).subs,{':'}) % obj.high{n}(:)
                % Extract all the values from a decomposition level
                out = [];
                for m = 1:length(obj.high{S(2).subs{1}})
                    out = cat(1, out, obj.high{S(2).subs{1}}{m}(:));
                end
            else
                out = builtin('subsref',obj,S);
            end
        end
        
        function obj = min(alpha,obj)
            if ~isa(alpha, 'waveletObj') && isa(obj, 'waveletObj')
                % case scalar * wavelet
                obj.low = min(alpha, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = min(alpha, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = min(alpha, obj.extra);
                end
            elseif isa(alpha, 'waveletObj') && ~isa(obj, 'waveletObj') 
                % case wavelet * scalar
                obj = min(obj,alpha);
            elseif isa(alpha, 'waveletObj') && isa(obj, 'waveletObj') 
                % case wavelet * wavelet
                obj.low = min(alpha.low, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = min(alpha.high{i}{m}, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = min(alpha.extra, obj.extra);
                end
            end
        end
        
        function obj = max(alpha,obj)
            if ~isa(alpha, 'waveletObj') && isa(obj, 'waveletObj')
                % case scalar * wavelet
                obj.low = max(alpha, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = max(alpha, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = max(alpha, obj.extra);
                end
            elseif isa(alpha, 'waveletObj') && ~isa(obj, 'waveletObj') 
                % case wavelet * scalar
                obj = max(obj,alpha);
            elseif isa(alpha, 'waveletObj') && isa(obj, 'waveletObj') 
                % case wavelet * wavelet
                obj.low = max(alpha.low, obj.low);
                for i = 1:length(obj.high)
                    for m = 1:length(obj.high{i})
                        obj.high{i}{m} = max(alpha.high{i}{m}, obj.high{i}{m});
                    end
                end
                if ~isempty(obj.extra)
                    obj.extra = max(alpha.extra, obj.extra);
                end
            end
        end
        
        function wavplot = imagesc(obj)
            % Determine if decimated
            obj2 = waveletObj(gather(obj));
            isDec = any(size(obj2.low) ~= size(obj2.high{1}{1}));
            obj2 = abs(obj2);
            switch ndims(obj2.low)
                case 1
                    error('TODO: define 1D wavelet plot');
                case 2
                    if isDec
                        wavplot = obj2.low;
                        % Determine levels
                        J = [1 1];
                        for j=2:length(obj2.high)
                            J = J + (size(obj2.high{j}{1}) < size(obj2.high{j-1}{1}));
                        end
                        if length(obj2.high{2})<3
                            % Wavelet level for one dim must be 0
                            [~,minInd] = min(J);
                            J(minInd) = 0;
                        end
                        for j=length(obj2.high):-1:1
                            % Loop through levels
                            if J(1) >= j
                                wavplot = cat(1,wavplot,obj2.high{j}{1});
                            end
                            sub = [];
                            if all(J>=j)
                                sub = cat(1,obj2.high{j}{2},obj2.high{j}{3});
                            elseif J(2) >= j
                                sub = obj2.high{j}{1};
                            end
                            wavplot = cat(2,wavplot,sub);
                        end
                    else
                        ml = 1;
                        for j=length(obj2.high):-1:1
                            ml = max(ml,length(obj2.high{j}));
                        end
                        wavplot = cat(2,obj2.low,zeros([size(obj2.low,1) (ml-1)*size(obj2.low,2)]));
                        for j=length(obj2.high):-1:1
                            sub = cell2mat(obj2.high{j});
                            if length(obj2.high{j}) < ml
                                sub = cat(2,sub, zeros([size(obj2.low,1) (ml-length(obj2.high{j}))*size(obj2.low,2)]));
                            end
                            wavplot = cat(1,wavplot,sub);
                        end
                    end
                    wavplot = abs(wavplot);
                    imagesc(wavplot);
                    colormap('gray'); axis off; axis image;
                    colorbar;
                case 3
                    % Plot slices in two orthogonal dimensions. Use subplot per level
                    wavplot = cat(1, obj2.low(:,:,round(end/2)+1), permute(squeeze(obj2.low(round(end/2)+1,:,:)), [2 1]) );
                    subplot(length(obj2.high)+1,1,1)
                    wavplot = abs(wavplot);
                    imagesc(gather(wavplot)), colormap('gray')
                    axis off; axis image;
                    
                    np = 2;
                    for j=length(obj2.high):-1:1
                        wavplot = cell2mat(obj2.high{j});
                        wavplot = cat(1, wavplot(:,:,round(end/2)+1), permute(squeeze(wavplot(round(end/2)+1,:,:)), [2 1]));
                        wavplot = abs(wavplot);
                        subplot(length(obj2.high)+1,1,np)
                        imagesc(gather(wavplot)), colormap('gray')
                        axis off; axis image;
                        np = np+1;
                    end
                otherwise
                    error('TODO: define higher dimensional wavelet plot');
            end
        end
        
	end
end

