function out = polysort(in)
% Sorts indices of a polynomial ROI so that they go in order around the ROI
% 
% (c) Corey Baron 2010

rest = in(2:end,:);
out = in;
out(1,:) = in(1,:);
for n=1:size(in,1)-1
    lengths = sqrt(sum((rest - repmat(out(n,:),size(rest,1),1)).^2,2));
    [~,inds] = sort(lengths);
    out(n+1,:) = rest(inds(1),:);
    if n < size(in,1)-1
        rest = rest(1:length(rest) ~= inds(1),:);
    end
end