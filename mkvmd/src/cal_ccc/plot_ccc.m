clear;
clc;

filename = ['ccc.dat'];
data=load(filename);

clf;figure(1);
plot(data(:,1),data(:,2),'-o', 'linewidth',2);
 
legend('CCC','fontsize',20, 'location','southeast');
title('', 'fontsize',20);
xlabel('ns', 'fontsize',20);
ylabel('CCC', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:50:250];
% ax.YTick=[0:2:7];
 
print('PLOT_ccc','-dpng','-r500');