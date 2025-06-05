function [flow,cutside] = mincut(sourcesink,remain,typ)


if ~exist('typ','var');typ='Int';end


%   MINCUT   Graph Max-Flow/Min-Cut.
%   [FLOW,CUTSIDE] = MINCUT(SOURCESINK,REMAIN) returns a graph's maximum flow and minimum cut.
%   The graph is defined by two matrices given as the two function arguments:
%   
%   SOURCESINK: nx3 matrix where the first column has a node number and the other
%   two the correspondent weights of the arcs linking that node to the source and
%   sink respectively. There is one line for each node (n = total amount of nodes).
%   
%   REMAIN: this is an mx4 matrix containing the information about the links between
%   nodes and respective weights. Column 1 goes for starting node, column 2 for final 
%   node and the final two columns have the direct and inverse link weights respectively.     
%
%   The outputs are given by:
%   
%   FLOW: a scalar.
%   
%   CUTSIDE: nx2 matrix where the first column has a node number and the second column
%   has one of the following values: 
%
%   0 - the node stays at the source side of the CUT.
%   1 - the node stays at the sink side of the CUT.


%   This m-file calls a mex file that contains a software library that is a modification 
%   of the maxflow algorithm described in:
% 
% 	An Experimental Comparison of Min-Cut/Max-Flow Algorithms
% 	for Energy Minimization in Computer Vision.
% 	Yuri Boykov and Vladimir Kolmogorov.
% 	In Third International Workshop on Energy Minimization
% 	Methods in Computer Vision and Pattern Recognition, September 2001


if strcmp(typ,'Int')
    sourcesinkg = int32(sourcesink);
    remaing = int32(remain);
elseif strcmp(typ,'Flo')
    sourcesinkg = single(sourcesink);
    remaing = single(remain);
else
    sourcesinkg = double(sourcesink);
    remaing = double(remain);
end
sourcesinkg(:,1)=sourcesinkg(:,1)-1;
remaing(:,1:2)=remaing(:,1:2)-1;
    
[flowg,cutsideg] = mf2(sourcesinkg,remaing);

flow = double(flowg); 
cutside = double(cutsideg);

% End