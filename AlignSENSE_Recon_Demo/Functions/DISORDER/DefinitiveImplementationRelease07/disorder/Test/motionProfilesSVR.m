T=svr.TVHist;
for m=4:4%size(T,6)
    figure
    for n=1:size(T,4)
        plot(180*T(:,1,1,n,m)/pi)
        hold on
    end
end

for v=1:svr.NV
    T=svr.TPHist{v};
    for m=4:4%size(T,6)
        figure
        for n=1:size(T,5)
            plot(180*T(:,1,1,1,n,m)/pi)
            hold on
        end
    end
end