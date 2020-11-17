clear;
clc;

PREFIX={'rmsd','rgyr'};

for jj=1:2
    ff=figure(jj);clf;
    cc=getPyPlot_cMap('Set1',9);
    for ii=1:3
        data=load(['res_',PREFIX{jj},'/',num2str(ii),'.dat']);
        hold on;
        plot(data(:,1)/10,data(:,2), 'linewidth',2,'color',cc(ii,:));
    end


    legend({'Trial 1' 'Trial 2' 'Trial 3'}, 'fontsize',28,'location','best');
    title('', 'fontsize',28);
    xlabel('', 'fontsize',28);
    ylabel('', 'fontsize',28);

    ax=gca;
    ax.FontSize=28;
    % xlim([0 100]);
    % ax.XTick=[0:50:100];
    % ylim([0 100]);
    % ax.YTick=[0:50:100];

%     print(['PLOT_',PREFIX{jj}],'-dpng');
end