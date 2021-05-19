clear;
clc;

filename = {'.dat' '.dat'};
f=figure(1);clf;
f.Position = [0 0 1024 1024];
% cc=getPyPlot_cMap('Set1',2);
for ii=1:2
    data=load(filename{ii});
    hold on;plot(data(:,1),data(:,2),'linewidth',2);
end

standard_axis() 
 
% print('PLOT_','-dpng');

function standard_axis() 
    ax=gca;
    ax.FontSize=48;
    % legend({'' ''},'fontsize',28,'location','best');
    % title('','fontsize',28);
    % xlabel('','fontsize',28);
    % ylabel('','fontsize',28);
    % xlim([0 100]);
    % ax.XTick=[0:50:100];
    % ylim([0 100]);
    % ax.YTick=[0:50:100];
end