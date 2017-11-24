clc;
clear;

filename={'avg.dat'};

data=load(filename{1});
figure(1);clf;
bb=bar(data(:,1),data(:,2),.7);
set(bb,'facecolor',lines(1));
set(bb,'edgecolor',lines(1));
set(bb,'linewidth',1);

set(gca, 'LineWidth',2);
ax=gca;
ax.FontSize=48;
xlim([251 361]);
% ax.XTick=[0:50:250];
% ylim();
ax.YTick=[0:.2:1];

fig = gcf;
fig.PaperPositionMode = 'auto';
str_pic=['PLOT_bar'];
print(str_pic,'-dpng','-r0');