close all;
clear;
clc;
name={'zaxis'};

kk=0;

% Time
% either set ns_start or frame_start
psperframe=10;
f2ns=psperframe/1000;
ns_eq=60;
frame_eq=ns_eq/f2ns;

for nn=1:2
    figure(nn);clf;
    data=[];
    data_all=[];
    for window=[583 595 607 619 631 643 655 667 679 691 703 715 727 739 751]
        filename = ['us-z',num2str(window),'/output_z/sel',num2str(nn),'.dat'];
        data_raw=load(filename);

        data=[data data_raw(end-frame_eq:end,2)];
        data_hist=histcounts(data_raw,[0:1:180],'Normalization', 'probability');
%         data_hist=histcounts(data_raw,[0:1:180]);
        
        data_all=[data_all;data_hist];
    end
%     histogram2([5.83 5.95 6.07 6.19 6.31 6.43 6.55 6.67 6.79 6.91 7.03 7.15 7.27 7.39 7.51],data);
    bb=bar3([58.3 59.5 60.7 61.9 63.1 64.3 65.5 66.7 67.9 69.1 70.3 71.5 72.7 73.9 75.1],data_all,1);
% %     set(bb,'facecolor',lines(1));
% %     set(bb,'edgecolor',lines(1));
    set(bb,'linestyle','none');
%     set(bb,'linewidth',1);
    colormap jet;
    cc=colorbar;
    cc.Ticks=[];
    for k = 1:length(bb)
        zdata = bb(k).ZData;
        bb(k).CData = zdata;
        bb(k).FaceColor = 'interp';
    end
    view(2);
    grid on;
    ax=gca;
%     ax.Color='w';
%     ax.XColor='w';
%     ax.YColor='w';
%     ax.TickLength=[.1 .5];
    ax.LineWidth=2;
    ax.FontSize=36;
    set(gca,'ticklength',3*get(gca,'ticklength'));
%     set(gca,'tickdir','out');
    ax.XTick=[0:20:180];
    ylim([58.3 75.1]);
    ax.YTick=[60 65 70 75];
    
    pbaspect([5 1 1]);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 24 9];
%     fig.PaperPositionMode = 'auto';
    set(gca,'LooseInset',get(gca,'TightInset'))
    str_pic=['PLOT_xaxis_pocket',num2str(nn)];
    print(str_pic,'-dpng','-r0');
end

%         legend({'PH1'},'fontsize',20, 'location','best');
%         title(['Pocket ',num2str(nn),' with Z'], 'fontsize',20);
% xlabel('ns', 'fontsize',20);
%         ylabel(str_y{jj}, 'fontsize',20);



