function y = simple_raw_filter(x,dim,startIdx,endIdx,typ,alpha,display)
%simple_raw_filter applies linear or exponential attenuation in a specific
%region of kspace
%   * x is the to-be-filtered kspace data
%   * dim is the dimention index of applying the attenuation filter
%   * startIdx is the start index of the line of the attenuated region in
%   {dim}
%   * endIdx is the end index of the line of the attenuated region in {dim}
%   * {typ} is the attenuation filter type 'linear' or 'exp' (for
%   expoential)
%   ** alpha is the exponential parameter for exponential filter
%   ** {display} 1- display the plot for appiled filter
%
%   Zihan Ning 20-01-2025

if nargin < 6; alpha=[];end
if nargin < 7; display = 0; end 
if contains(typ, 'exp') && isempty(alpha); alpha = 5; end % default value for exponential model

numSlices = endIdx - startIdx + 1;

% generate attenuation model
if contains(typ, 'exp')
    temp_filter = exp(-alpha * (0:numSlices-1)/(numSlices-1));
else % linear in defualt
    temp_filter = linspace(1, 0, numSlices); % [1, ..., 0]
end

% match the demension
if length(size(x)) < dim
    fprintf('Entered dimension to apply raw filter is larger than the whole dimension!!!! \n');
    fprintf('Raw filter is not applied \n');
    y = x;
    return;
else % be able to apply the filter
    apply_dim = ones(1,length(size(x)));
    apply_dim(dim) = numSlices;
end
atten_filter = reshape(temp_filter, apply_dim);

% display the filter & signal before annuation
if display == 1
    figure; 
    plot(startIdx:endIdx, temp_filter, '-o');
    xlabel('3rd Dimension Index');
    ylabel('Attenuation Factor');
    title('Linear Attenuation Across Specified Region');
    grid on;
end

% apply the attenuation filter
y = x; 
numDims = ndims(x);
idx = repmat({':'}, 1, numDims);
idx{dim} = startIdx:endIdx;
y(idx{:}) = y(idx{:}) .* atten_filter;
fprintf('Attenuation %s filter applied to %s dimension from %s to %s idx. \n',typ,num2str(dim),num2str(startIdx),num2str(endIdx));

