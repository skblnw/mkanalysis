clear;
clc;

bin=400;
window=15;

figure(1);clf;
filename = ['output_pmf/histo-50ns-400.xvg'];
fileID = fopen(filename);
formatSpec='';
for ii=1:window+1
    formatSpec = [formatSpec '%f'];
end
data = textscan(fileID, formatSpec, 'CommentStyle', {'#','@'}, 'HeaderLines', 15);
fclose(fileID);

for ii=2:window+1
    hold on;
    total=sum(data{ii});
    plot(data{1},data{ii}./total,'linewidth',2);
%     plot(data{1},data{ii},'linewidth',2);
end

% legend([],'fontsize',24,'location','southeast');
title('Histogram','fontsize',20);
xlabel('Separation (nm)','fontsize',20);
% ylabel('PMF (kcal/mol)','fontsize',20);
ax=gca;
ax.FontSize=20;
xlim([5.6 7.9]);
ax.XTick=[5.7:0.3:7.8];

% print(['FIG_hist_50ns_normalized'],'-dpng');