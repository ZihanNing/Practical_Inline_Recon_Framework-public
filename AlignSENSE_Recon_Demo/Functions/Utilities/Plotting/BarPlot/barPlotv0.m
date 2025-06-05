
function [] = barPlot(error, n_max, scans,Ylab)


%%% set variables
co=[0.8500    0.3250    0.0980;         
    0         0.4470    0.7410;    
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840
    0.4470    0.7410    0
    0.7410    0         0.4470
    0         0         0];
FontSize = 16;

title_text = '';
legendNames = {}; 

numScans = size(error,1); 
numSH = size(error, 2);
if size(error,3)<2;error = cat(3, error , zeros(size(error))) ; end

indx = 0:800; % Large enough if numbers of levels increases

figure('Name','','color', 'w')
p= {};
list_plots = zeros (1, numSH);

for i =1:numSH
        indexes =  indx(1:numScans)*numSH +  indx(1:numScans)+i ;  % Voor 6 segments

        mean = error (:,i,1);
        stdv = error(:,i,2); 
        
        p{i} = bar(indexes,mean,1/(numSH+1), 'FaceColor',co(i,:));hold on
        list_plots(i) = p{i};
        
        errorbar(indexes, mean,stdv,'LineStyle','none','Color',[0 0 0]);

        if i==1; legendNames = cat(2, legendNames, {('$Linear\ only$')});else legendNames = cat(2, legendNames, {sprintf('$n_{max} = %d$',n_max(i-1)) });end
end

hold off
indexes_axis = indx(1:numScans)*numSH + indx(1:numScans)+(numSH+1)/2 ;
for oo = 1:numScans; labels_axis{oo} =  num2str(scans(oo));end

set(gca,'XTick',indexes_axis,'XTickLabel',labels_axis);
set(gca, 'FontSize',FontSize, 'TickLabelInterpreter','latex');

title( title_text, 'interpreter','latex','FontSize', FontSize)
rangeLow = error(:,:,1)-error(:,:,2);
rangeUp = error(:,:,1)+error(:,:,2);
range = gather(1.1* [ min(min(rangeLow(:)),0)  max(rangeUp(:)) ]    );
xlim([0 numScans*numSH+1+numScans-1 ])
ylim([range(1) range(2)])  

AX=legend([list_plots],legendNames,'Location','NorthEast','FontSize',12,'interpreter','latex');    
LEG = findobj(AX,'type','text');
set(AX,'Position',get(AX,'Position')+[0 0 0 0])%To the right and to the top         
                
ylabel(Ylab, 'interpreter', 'latex', 'FontSize', FontSize)
xlabel('Scan', 'interpreter', 'latex', 'FontSize', FontSize)
