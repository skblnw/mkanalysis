clear;
clc;


filename = ['rmsd_npt.dat'];
data=load(filename);
psperframe=100;
data(:,1)=data(:,1)*psperframe/1000;

eq_ns=100;
eq_frame=eq_ns*1000/psperframe;

str_tab=['TABLE_avg.txt'];
id_outfile = fopen(str_tab,'w');

figure(1);clf;
plot(data(:,1),data(:,2),'linewidth',2);
str_leg=['RMSD'];
fprintf(id_outfile,'%s\t%.2f (%.2f)\n',str_leg,mean(data(end-eq_frame:end,2)),std(data(end-eq_frame:end,2)));

legend(str_leg,'fontsize',20, 'location','best');
title('', 'fontsize',20);
xlabel('ns', 'fontsize',20);
ylabel('RMSD (A)', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% xlim([0 1200]);
% ax.XTick=[0 500 1000];
% ylim([0 3]);
% ax.YTick=[0:1:3];
 
% print('PLOT_rmsd','-dpng','-r500');

fclose(id_outfile);