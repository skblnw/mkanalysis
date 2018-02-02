clear;
clc;

name={'FWCSS', 'AWCSS', 'BWCSS'};

for jj=1:3

    figure(jj);clf;
    wcss_final=[];
    for ii=1:372
        filename = [name{jj},'_cg4_init',num2str(ii),'_iter.dat'];
        data=load(filename);
        wcss_final=[wcss_final data(end)];

        nn=size(data,1);
        hold all;
        if nn == 1
            pp=plot(0,data,'-*','linewidth',1.5,'markersize',6,'MarkerFaceColor','b');
            pp.Color(4)=0.2;
        else
            pp=plot(linspace(0,1,nn),data,'-*','linewidth',1.5,'markersize',6,'MarkerFaceColor','b');
            pp.Color(4)=0.2;
        end
    end

    [M, I] = min(wcss_final);
    fprintf('%s: the smallest is %d\n', name{jj}, I);

    for ii=279
        filename = [name{jj},'_cg4_init',num2str(ii),'_iter.dat'];
        data=load(filename);
        nn=size(data,1);
        hold all;
        if nn == 1
            plot(0,data,'-*k','linewidth',3,'markersize',10,'MarkerFaceColor','b');
        else
            plot(linspace(0,1,nn),data,'-*k','linewidth',3,'markersize',10,'MarkerFaceColor','b');
        end
    end

    % legend('','fontsize',20, 'location','best');
    title('', 'fontsize',20);
    xlabel('Steps', 'fontsize',20);
    ylabel('WCSS', 'fontsize',20);

    grid on;
    ax=gca;
    ax.FontSize=20;
    xlim([-.1 1.1]);
    ax.XTick=[];
    % ylim([26.5 29.5]);
    % ax.YTick=[0:5:100];

    print(['PLOT_profile_',name{jj}],'-dpng');

end