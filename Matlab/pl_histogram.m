clear;
clc;

PREFIX={'hbonds'};

data={};
for jj=1:1
    figure(jj);clf;
    sp1=subplot(1,2,1);
    cc=getPyPlot_cMap('Set1',9);
    for ii=1:3
        data{ii}=load(['res_',PREFIX{jj},'/',num2str(ii),'.dat']);
        hold on;
        plot(data{ii}(:,1)/10,data{ii}(:,24), 'linewidth',1,'color',cc(ii,:));
    end


    legend({'Trial 1' 'Trial 2' 'Trial 3'}, 'fontsize',28,'location','best');
    title('', 'fontsize',28);
    xlabel('', 'fontsize',28);
    ylabel('', 'fontsize',28);

    ax=gca;
    ax.FontSize=28;
    % xlim([0 100]);
    % ax.XTick=[0:50:100];
    ylim([10 35]);
    % ax.YTick=[0:50:100];
    
    sp2=subplot(1,2,2);
    datahist=[];
    for ii=1:2
        datahist=[datahist;data{ii}(size(data{ii},1)/2:end,24)];
    end
    histogram(datahist,'Normalization','pdf','facecolor','k');
    ax=gca;
    xlim([10 35]);
%     ylim([0 .1]);
    ax.XTick=[];
    ax.YTick=[];
    set(gca,'xdir','reverse')
    camroll(-90);
    
    ah1=set(sp1,'Units','normalized','Position',[.2 .1 .5 .8]);
    ah2=set(sp2,'Units','normalized','Position',[.75 .1 .2 .8]);

%     print(['PLOT_',PREFIX{jj}],'-dpng');
end