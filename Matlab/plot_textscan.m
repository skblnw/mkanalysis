clear;
clc;

filename = [''];
% data=load(filename);
% filename = ['',num2str(i), ''];
fileID = fopen(filename);
formatSpec = '%s %f';
data = textscan(fileID, formatSpec, 'CommentStyle', {'#' ,'@' }, 'HeaderLines', 0);
% For excluding data
% frame_start=850;
% data(frame_start:end,:)=[];

figure(1);clf;
plot(data(:,1),data(:,2), 'linewidth',2);
 
legend('','fontsize',20, 'location','best');
title('', 'fontsize',20);
xlabel('ns', 'fontsize',20);
ylabel('', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% xlim([0 1000]);
% ax.XTick=[0:50:250];
% ylim([0 10]);
% ax.YTick=[0:2:7];
 
% print('PLOT_rmsd','-dpng');