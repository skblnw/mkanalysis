clc;
clear;

filename={'reslist_complex'};

data=load(filename{1});
figure(1);clf;
bb=bar(data(:,1)/1000,.7);
set(bb,'facecolor',lines(1));
set(bb,'edgecolor',lines(1));
set(bb,'linewidth',1);

set(gca, 'LineWidth',2);
ax=gca;
ax.FontSize=24;
xlabel("Residue");
ylabel("Counts");
xlim([600 1000]);
% ax.XTick=[0:50:250];
% ylim();
% ax.YTick=[0:.2:1];

fig = gcf;
fig.PaperPositionMode = 'auto';
str_pic=['PLOT_bar_750_1000'];
% print(str_pic,'-dpng','-r0');