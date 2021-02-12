clear;
clc;

filename = {'.dat' '.dat'};
figure(1);clf;
% cc=getPyPlot_cMap('Set1',2);
for ii=1:2
    data=load(filename{ii});
    hold on;
    plot(data(:,1),data(:,2), 'linewidth',2);
end
 
legend({'' ''}, 'fontsize',28, 'location','best');
title('', 'fontsize',28);
xlabel('', 'fontsize',28);
ylabel('', 'fontsize',28);

ax=gca;
ax.FontSize=28;
% xlim([0 100]);
% ax.XTick=[0:50:100];
% ylim([0 100]);
% ax.YTick=[0:50:100];
 
% print('PLOT_','-dpng');

