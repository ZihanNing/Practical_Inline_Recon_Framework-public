function X=bldiag2arrGPU(C,indBl)

%BLDIAG2ARRGPU   Converts (forwards and backwards) a cell-based block diagonal
%representation to a lexicographically ordered array representation. The direction 
%of conversion is automatically determined by checking the type of C
%   X=BLDIAG2ARRGPU(C,INDBL)
%   * C are the set of blocks arranged in cells / lexicographically ordered 
%   array
%   * INDBL are the sizes of the blocks
%   * C is the lexicographically ordered array / set of blocks arranged in cells
%

di=iscell(C);%Direction of conversion

NC=length(indBl);
if di            
    for n=1:NC
        NCA=size(C{n});NCA(end+1:3)=1;
        C{n}=reshape(C{n},[prod(NCA(1:2)) 1 NCA(3:end)]);
    end
    X=cat(1,C{:});
else
    X=cell(1,NC);
    indBEnd=cumsum(indBl.^2);
    indBSta=[1 indBEnd(1:end-1)+1];
    
    for n=1:NC
        X{n}=dynInd(C,indBSta(n):indBEnd(n),1);
        NCA=size(X{n});NCA(end+1:3)=1;
        NCA(1:2)=gather(indBl(n));
        X{n}=reshape(X{n},NCA);
    end
end
