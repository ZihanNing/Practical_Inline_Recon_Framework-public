
sour='/home/lcg13/Data/rawDestin/ReconstructionsDebug06';
path='2018_03_26/DI_4422';
fil=ls(sprintf('%s/%s/An-S2/*SVR.mat',sour,path));

load(fil(1:end-1));

%return

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

for v=1:svr.NV
    for m=4:4
        figure
        for p=1:svr.P(v)
            T=svr.TEHist{v}{p};
            for n=1:size(T,6)
                plot(180*T(:,1,1,1,1,n,m)/pi)
                hold on
            end
        end
    end
end
            
        