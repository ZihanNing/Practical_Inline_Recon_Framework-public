load IndIntNew
codedSize=IndIntNew.coded_size;
offset=IndIntNew.offset;

NP=length(codedSize);
codedSize=codedSize([end 1:end-1]);
fixOf=find(codedSize~=0);%Fixed offset
codedSizeFilled=codedSize;
offsetModified=IndIntNew.offset;
NFO=length(fixOf);

for n=1:NFO
    if n~=NFO
        indZ=fixOf(n)+1:fixOf(n+1)-1;
    else
        indZ=fixOf(n)+1:NP;
    end
    if ~isempty(indZ)
        codedSizeFilled(indZ)=codedSize(fixOf(n));
        offsetModified(indZ)=offset(indZ)+cumsum(uint64(codedSizeFilled(indZ)));
        if n~=NFO
            codedSizeFilled(fixOf(n+1))=codedSizeFilled(fixOf(n+1))-sum(uint64(codedSizeFilled(indZ)));            
        else
            codedSizeFilled(1)=codedSizeFilled(1)-sum(uint64(codedSizeFilled(indZ)));
        end
        %if length(indZ>1) && n~=1
        %    codedSizeFilled(fixOf(n))=offsetModified(fixOf(n))-offsetModified(fixOf(n)-1);
        %end
    end
end
codedSizeFilled=codedSizeFilled([2:end 1]);
% figure
% plot(codedSize(1:32))
% figure
% plot(codedSizeFilled(1:32))
% figure
% plot(offset(1:32))
% figure
% plot(offsetModified(1:32));

%figure
%plot(codedSize(end-31:end))

figure
plot(codedSizeFilled(end-31:end))

%figure
%plot(IndIntNew.size(end-31:end))

figure
plot(codedSizeFilled(1:32))

%figure
%plot(IndIntNew.size(1:32))


%figure
%plot(offset(end-31:end))
%figure
%plot(offsetModified(end-31:end));
% 
% figure
% plot(codedSize)
% figure
% plot(codedSizeFilled)
% figure
% plot(offset)
% figure
% plot(offsetModified);
% 
%     