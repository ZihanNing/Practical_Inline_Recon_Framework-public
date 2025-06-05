function visArray(x)

%VISARRAY   Visualizes a multidimensional array (up to 4D) by concatenating
%the third dimension in the columns and the fourth dimension in the rows
%   VISARRAY(X)
%   * X is the array to be visualized
%

if ~isempty(x)
    fprintf('Displaying array of size%s\n',sprintf(' %d',size(x)));
    x=squeeze(x);
    nd=numDims(x);
    assert(nd<=4,'Visualization of arrays is limited to 4 dimensions');
    N=size(x);N(end+1:4)=1;
    x=permute(x,[1 4 2 3]);
    x=reshape(x,[prod(N([1 4])) prod(N([2 3]))]);
    if ~any(imag(x(:)))
        figure
        imshow(x,[]);
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);     
    else
        figure
        imshow(abs(x),[])
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);     
        figure
        imshow(angle(x),[])
        set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);     
    end
end
 
