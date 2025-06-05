

function [AReshaped,pReshaped,idxReshape, ptemp] = reshapeP(A,p,type,idx,isOffset)

if nargin<2 || isempty(p);p=[];end
if nargin<3 || isempty(type);type='mat2vec';end
if nargin<4 || isempty(idx);if strcmp(type,'vec2mat');error('reshapeP:: Need indices for reshaping from vector to matrix');end;end
if nargin<5 || isempty(isOffset);isOffset=0;end%By default

%%% DEDUCE PARAMETER SIZE
if strcmp(type,'mat2vec')
   ndG = size(A,1);
   if ~isempty(p) && ~isempty(A);isOffset= size(A,2)~=size(p,1);end
   if ~isempty(p); NCha = size(p,1)-isOffset;else;NCha = size(A,2)-isOffset;end
elseif strcmp(type,'vec2mat')
   ndG = max(idx);
   NCha = length(idx)/ndG-isOffset;%isOffset must be provided in the input
else
    error('reshapeP:: Type not supported.')   
end

if ~isempty(p)
    if strcmp(type,'mat2vec')
       ndS = size(p,2);
    elseif strcmp(type,'vec2mat')
       ndS = size(p,3);
    end
end

%%% RESHAPE
%Calibration matrix
if strcmp(type,'mat2vec')
    AReshaped=[];idxReshape=[];
    for i=1:ndG
        AReshaped=cat(1,AReshaped,A(i,:).');
        idxReshape = cat(1,idxReshape, i*ones([length(A(i,:)) 1]));
    end
elseif strcmp(type,'vec2mat')
   AReshaped = zeros([ndG NCha+isOffset],'like',A);for i=1:ndG; AReshaped = dynInd(AReshaped,i,1,dynInd(gather(A),idx==i,1).');end 
   idxReshape=[];
end

%Pilot Tone signal
if  ~isempty(p)
    if strcmp(type,'mat2vec') 
        pReshaped=zeros([ndG ndG*(NCha+isOffset) ndS], 'like', p);
        if isOffset; p = padArrayND(p,[1 0],1,'post',1);end%Offset in calibration
        for i=1:ndG
            for ii=1:ndS
                pReshaped = dynInd(pReshaped,{i,idxReshape==i, ii},1:3,dynInd(p,{ii},2).');
            end 
        end
    elseif strcmp(type,'vec2mat') 
        pReshaped=zeros([NCha ndS], 'like', p);
        for i=1:ndS
            pReshaped = dynInd(pReshaped, i, 2, dynInd(p, {1,1:NCha,i},1:3));
        end
    end
else
    pReshaped=[];
end

%------------- temp
if ~isempty(pReshaped) &&  strcmp(type,'mat2vec') 
    assert(ndG==size(pReshaped,1),'temp error')
    ptemp =zeros([size(pReshaped,3)*ndG size(pReshaped,2)],'like',p);
    for i=1:size(pReshaped,3); ptemp= dynInd(ptemp, (1+(i-1)*ndG):i*ndG,1, dynInd(pReshaped,i,3));end
else
    ptemp=[];
end
    
    
    
