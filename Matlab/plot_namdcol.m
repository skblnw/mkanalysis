clear;
clc;
close all;

name={'colvars' 'ener' 'temp' 'box'};

cc=lines(4);
for window=[583 595 607 619 631 643 655 667 679 691 703 715 727 739 751 763]
    for ii=1:4
        filename = ['us-z',num2str(window),'-',name{ii},'.dat'];
        data=load(filename);

        figure(ii);clf;
        plot(data(:,1),data(:,2), 'linewidth',2,'color',cc(ii,:));

%         legend('','fontsize',20, 'location','best');
        title(['us-z',num2str(window),'-',name{ii}], 'fontsize',20);
        xlabel('ns', 'fontsize',20);
        ylabel('', 'fontsize',20);

        ax=gca;
        ax.FontSize=20;
        % xlim([0 1000]);
        ax.XTick=[0:20:500];
        % ylim([0 10]);
        % ax.YTick=[0:2:7];

        print(['PLOT_us-z',num2str(window),'-',name{ii}],'-dpng');
    end
end