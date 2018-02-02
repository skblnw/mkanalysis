clear;
clc;

figure(1);clf;

filename = ['bm-result.dat'];
data=load(filename);
data=dlmread(filename, ' ', 1, 0);

cpu_per_node=24;

hold on;
plot(data(:,1)*cpu_per_node,data(:,2),'-o','linewidth',2,'markersize',10);

grid on;
grid minor;
ax=gca();
ax.XGrid='off';
ax.XMinorGrid='off';
ax.LineWidth=1;
legend(['Benchmark (',cpu_per_node,' cpu/node'], ...
    'fontsize',8,'location','northwest');
title('', 'fontsize',20);
xlabel('Number of cores', 'fontsize',20);
ylabel('ns/day', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:50:250];
% ax.YTick=[0:2:7];
 
% print('benchmark-compare-all','-dpng','-r500');