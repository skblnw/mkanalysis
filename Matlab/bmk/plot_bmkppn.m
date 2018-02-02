clear;
clc;

figure(1);clf;

filename = ['combo-12.dat'];
data=load(filename);
data(:,2)=data(:,2)./data(end,2)*100;
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);

filename = ['tj-1.dat'];
data=load(filename);
data(:,2)=data(:,2)./data(end,2)*100;
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);

filename = ['tj-4.dat'];
data=load(filename);
data(:,2)=data(:,2)./data(end,2)*100;
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);

filename = ['gz-actin-1.dat'];
data=load(filename);
data(:,2)=data(:,2)./data(end,2)*100;
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);

filename = ['gz-water-1.dat'];
data=load(filename);
data(:,2)=data(:,2)./data(end,2)*100;
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);

legend({'Combo (K20) 12 nodes', ...
        'Tianjin (M2050) 1 node', ...
        'Tianjin (M2050) 4 nodes', ...
        'Guangzhou (M2050) 500k atoms 1 node', ...
        'Guangzhou (M2050) 50k atoms 1 node'}, ...
    'fontsize',14, 'location','southeast');
title(sprintf('NAMD2.11\n# of cores per node for GPU'), 'fontsize',20);
xlabel('Number of cores', 'fontsize',20);
ylabel('Benchmark (%)', 'fontsize',20);

% xlim([1 13]);
% ylim([0 120]);
grid on;
grid minor;
ax=gca;
ax.XGrid='off';
ax.XMinorGrid='off';
ax.LineWidth=1;
% ax.MinorGridLineStyle='--';
ax.FontSize=20;
% ax.XTick=[8 12];
% ax.YTick=[0:20:100];
 
print('benchmark-gpu-ppn','-dpng','-r500');