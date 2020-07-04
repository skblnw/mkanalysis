clear;
clc;

filename = {'vmi534-3880_ups1-mdm35.dat' 'vmi532-3881_ups1-mdm35_h33hsp.dat'};
figure(1);clf;hold all;
data1=load(filename{1});
% plot(data(:,1)/150,data(:,2), 'linewidth',2);
data2=load(filename{2});
% plot(data(:,1),data(:,2), 'linewidth',2);

% edges=[0:.2:15];
% h=histcounts(data,edges,'Normalization', 'probability');
% plot(edges,[0 h]);

h1=histogram(data1(:,2),'Normalization','pdf');
hold on;
h2=histogram(data2(:,2),'Normalization','pdf');

h1.BinLimits=[min(min([h1.BinLimits h2.BinLimits])) max(max([h1.BinLimits h2.BinLimits]))];
h2.BinLimits=h1.BinLimits;
h1.NumBins=min([h1.NumBins h2.NumBins]);
h2.NumBins=h1.NumBins;

legend({'' ''}, 'fontsize',28,'location','best');
title('', 'fontsize',28);
xlabel('', 'fontsize',28);
ylabel('Probability', 'fontsize',28);

ax=gca;
ax.FontSize=28;
% xlim([0 100]);
% ax.XTick=[0:50:100];
% ylim([0 100]);
ax.YTick=[];
 
% print('PLOT_','-dpng');
