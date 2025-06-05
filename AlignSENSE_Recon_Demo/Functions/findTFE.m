function [TFE] = findTFE(TW,dim)

y = TW.image(:,1,:,:,1,1,1,1,1,1);
y = squeeze(y);

y = dynInd(y,centerIdx(size(y,1)),1);

if dim==2
    y = sum(abs(y),3);
elseif dim==3
    y = sum(abs(y),2);
end
    
    
y = squeeze(y);
TFE = sum( single(y>0));