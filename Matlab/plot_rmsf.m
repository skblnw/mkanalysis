clear;
clc;

figure(1);clf;

filename = ['rmsf_protein.dat'];
data=load(filename);
hold on;
plot(data(:,1),data(:,2), 'linewidth',2);
 
legend({'protein'},'fontsize',20, 'location','best');
% title('Ups1', 'fontsize',20);
xlabel('Residue', 'fontsize',20);
ylabel('RMSF(A)', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:20:180];
% ax.YTick=[0:2:7];
 
% print('PLOT_rmsf','-dpng','-r500');