clear;
clc;

filename = ['rmsd_npt.dat'];
data=load(filename);
psperframe=100;
data(:,1)=data(:,1)*psperframe/1000;

figure(1);clf;
plot(data(:,1),data(:,2), 'linewidth',2);
 
legend('RMSD','fontsize',20, 'location','best');
title('', 'fontsize',20);
xlabel('ns', 'fontsize',20);
ylabel('RMSD (A)', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:50:250];
% ax.YTick=[0:2:7];
 
% print('PLOT_rmsd','-dpng','-r500');