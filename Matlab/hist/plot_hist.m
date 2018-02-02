clear;
clc;
% name1: name of data in one plot
name1={'cls1','cls2'};

% Histogram bin size
bin_size=.1;

% Y-axis label
str_y={'Angle (^{\circ})'};

str_tmp=['TABLE_mean.txt'];
id_outfile = fopen(str_tmp,'w');

figure(kk);clf;
% Counting data set
kk=0;
for nn=1:2
    kk=kk+1;
    filename = ['ang_',name1{nn},'.dat'];
    data=load(filename);
    
    % Plot hitogram
    % edges is unique for each plot
    edges=[120:bin_size:125];
    
    % set specific parameters for each data set
    % line color etc.
    switch nn
        case 1
            style=['-'];
            color=['b'];
            str_leg{kk}=[name{ii}];
        case 2
            style=['-'];
            color=['r'];
            str_leg{kk}=[name{ii}];
    end
    
    h=histcounts(data,edges,'Normalization', 'probability');
    hold on;
    pp(kk)=plot(edges,[0 h],style,'color',color,'linewidth',3,'markersize',10);
    
    % Write average
    fprintf(id_outfile,'%s\t%.2f (%.2f)\n',[name2{nn}],mean(data_hist),std(data_hist));

    % Write histogram data
    dlmwrite(['hist/DAT_hist-cls1-',name2{nn},'-mdff.dat'],data_hist);
end
fclose(id_outfile);

xlim([edges(1) edges(end)]);
%     ylim([0 0.2]);
ax=gca;
ax.FontSize=20;
ax.XTick=[edges(1) (edges(1)+edges(end))/2 edges(end)];

legend(p,str_leg,'fontsize',36,'location','northeast');

str_pic=['PLOT_hist'];
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10.24 9.6];
%     print(str_pic,'-dpng','-r0')