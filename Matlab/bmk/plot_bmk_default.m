clear;
clc;

filename = ['bm-result.dat'];
data=load(filename);
data=dlmread(filename);

cpu_per_node=24;

figure(1);clf;
hold on;
plot(data(:,1),data(:,2),'-bo','linewidth',2,'markersize',5);
ylabel('ns/day', 'fontsize',20);

figure(2);clf;
for row=1:size(data,1)
    data(row,3) = data(row,2)/data(row,1);
    data(row,4) = data(row,3)/data(1,2);
end
plot(data(:,1),data(:,4),'-ro','linewidth',2,'markersize',5);
ylabel('Efficiency', 'fontsize',20);

for ii=1:2
    figure(ii);
    grid on;
%     grid minor;
    ax=gca();
    ax.XGrid='off';
    ax.XMinorGrid='off';
    ax.LineWidth=2;
    %legend(['Benchmark (',cpu_per_node,' cpu/node'], 'fontsize',8,'location','northwest');
    title('~500k atoms, 16 cores/node', 'fontsize',20);
    xlabel('Number of nodes', 'fontsize',20);

    ax=gca;
    ax.FontSize=20;
    % ax.XTick=[0:50:250];
    % ax.YTick=[0:2:7];
end
 
print(figure(1), 'FIG_bmk_core16_omp4','-dpng','-r300');
print(figure(2), 'FIG_bmk_core16_omp4_strong','-dpng','-r300');