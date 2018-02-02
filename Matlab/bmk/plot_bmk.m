clear;
clc;

figure(1);clf;

filename = ['combo-new3.dat'];
data=load(filename);
hold on;
plot(data(:,1)*16,data(:,2),'-x','linewidth',2,'markersize',10);

if 0
    filename = ['titan-2.dat'];
    data=load(filename);
    hold on;
    plot(data(:,1)*16,data(:,2),'-x','linewidth',2);
else
    filename = ['titan-1.dat'];
    data=load(filename);
    hold on;
    plot(data(:,1)*16,data(:,2),'-x','linewidth',2);
end

filename = ['gz-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*24,data(:,2),'-x','linewidth',2);

filename = ['shenzhen-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1),data(:,2),'x','linewidth',2,'markersize',10);

filename = ['tianjin-gpu.dat'];
data=load(filename);
hold on;
plot(data(:,1),data(:,2),'-x','linewidth',2);

filename = ['combo-cpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*16,data(:,2),'-o','linewidth',2);

filename = ['college.dat'];
data=load(filename);
hold on;
plot(data(:,1)*16,data(:,2),'-o','linewidth',2);

filename = ['gz-cpu.dat'];
data=load(filename);
hold on;
plot(data(:,1)*24,data(:,2),'-o','linewidth',2);

filename = ['shenzhen.dat'];
data=load(filename);
hold on;
plot(data(:,1)*12,data(:,2),'-o','linewidth',2);

filename = ['tianjin-cpu.dat'];
data=load(filename);
hold on;
plot(data(:,1),data(:,2),'-o','linewidth',2);



grid on;
grid minor;
ax=gca();
ax.XGrid='off';
ax.XMinorGrid='off';
ax.LineWidth=1;
legend({'Combo (GPU-K20) 16-cores', ...
    'Titan (GPU-K20X AMD 6140 2.6GHz) 16-cores', ...
    'Guangzhou (GPU-M2050) 24-cores', ...
    'Shenzhen (GPU-C2050) 12-cores', ...
    'Tianjin (GPU-M2050) 12-cores', ...
    'Combo (CPU, E5-2660 2.2GHz) 16-cores', ...
    'College (CPU, E5-2670 2.60GHz) 16-cores', ...
    'Guangzhou (CPU, E5-2692 2.2GHz) 24-cores', ...
    'Shenzhen (CPU, X5650 2.66GHz) 12-cores', ...
    'Tianjin (CPU, X5670 2.93GHz) 12-cores'}, ...
    'fontsize',8, 'location','northwest');
title('', 'fontsize',20);
xlabel('Number of cores', 'fontsize',20);
ylabel('ns/day', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:50:250];
% ax.YTick=[0:2:7];
 
% print('benchmark-compare-all','-dpng','-r500');