clear;
clc;

filename = {'rmsf_ups1-pa_1.dat' 'rmsf_ups1-pa_2.dat'};
figure(1);clf;
for ii=1:2
    data=load(filename{ii});
    hold on;
    plot(data(:,1),data(:,2), 'linewidth',2);
end
 
legend({'Trial 1' 'Trial 2'},'fontsize',20, 'location','best');
title('Ups1-PA', 'fontsize',20);
xlabel('Residues', 'fontsize',20);
ylabel('RMSF (A)', 'fontsize',20);

ax=gca;
ax.FontSize=20;
xlim([0 169]);
% ax.XTick=[0:50:250];
ylim([0 3.5]);
% ax.YTick=[0:2:7];
 
print('PLOT_rmsf_ups1-pa','-dpng');