clear;
clc;

data=load('out');
figure(1);clf;
plot_corr_matrix(data);

xlabel('Residues','Fontsize',24);
ylabel('Residues','Fontsize',24);
% xticks([20:20:152]);
% 
% ax=gca();
% ax.FontSize=24;

colorbar;

print('PLOT_gc_r800-860','-dpng');