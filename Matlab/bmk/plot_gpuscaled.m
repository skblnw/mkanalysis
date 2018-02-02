clear;
clc;

figure(1);clf;

filename = ['combo-new2.dat'];
data=load(filename);
hold on;
plot(data(:,1)*12,data(:,2),'-x','linewidth',2);

filename = ['titan-1.dat'];
data=load(filename);
hold on;
plot(data(:,1)*12,data(:,2),'-x','linewidth',2);

filename = ['gz-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*8,data(:,2),'-x','linewidth',2);

filename = ['shenzhen-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*8/12,data(:,2),'x','linewidth',2,'markersize',10);

filename = ['tianjin-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*8/12,data(:,2),'-x','linewidth',2);

grid on;
grid minor;
ax=gca();
ax.XGrid='off';
ax.XMinorGrid='off';
ax.LineWidth=1;
legend({'Combo (GPU-K20) 12/16-cores', ...
    'Titan (GPU-K20X) 12/16-cores', ...
    'Guangzhou (GPU-M2050) 8/24-cores', ...
    'Shenzhen (GPU-C2050) 8/12-cores', ...
    'Tianjin (GPU-M2050) 8/12-cores'}, ...
    'fontsize',12, 'location','southeast');
title(sprintf('NAMD2.11\n516747 atoms\nGPU/CPU cores with max. performance'), 'fontsize',12);
xlabel('Number of cores', 'fontsize',20);
ylabel('ns/day', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:50:250];
% ax.YTick=[0:2:7];
 
print('benchmark-gpu-scaled','-dpng','-r500');