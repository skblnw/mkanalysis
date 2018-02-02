clear;
clc;
% name1: name of files
% name2: renaming name1 (if needed)
name1={'tetramer','ph-in-same-row'};
name2={'front','same'};

% Y-axis label
str_y={'Angle (^{\circ})'};

% Time (ns)
psperframe=10;
f2ns=psperframe/1000;
ns_start=10;
frame_start=ns_start/f2ns;

str_tmp=['TABLE_mean.txt'];
id_outfile = fopen(str_tmp,'w');
% Counting figure no.
kk=0;
for nn=1:2
    filename = ['ang_',name1{nn},'.dat'];
    data=load(filename);
    % modify data loaded (if needed)
    data=data(:,2:end);
    data_production=data(frame_start:end,:);

    % Re-arrange data_production to become data_hist
    data_hist=reshape(data_production,[],1);

    % Plot different columns first
    kk=kk+1;
    figure(kk);clf;
    for ii=1:length(data(1,:))
        hold on;
        pp=plot((1:length(data))*f2ns,data(:,ii),'--','linewidth',1);
        pp.Color(4)=0.4; % Transparency
        str_leg{ii}=['column ',num2str(ii)];
    end

    % Compute data_mean
    data_mean=mean(data,2);
    data_production_mean=data_mean(frame_start:end,:);
    % Re-arrange data_production_mean to become data_hist_mean
    data_hist_mean=reshape(data_production_mean,[],1);
	
    % Plot the average
    hold on;
    pp=plot((1:length(data_mean))*f2ns,data_mean(:,2),'-','linewidth',2);
    str_leg=[str_leg {'Average'}];
    legend(str_leg,'fontsize',24,'location','bestoutside');
    title('', 'fontsize',20);
    xlabel('ns', 'fontsize',20);
    ylabel('', 'fontsize',20);
    ax=gca;
    ax.FontSize=20;
%     str_pic=['hist/PLOT_cols-',name1{nn}];
%     print(str_pic,'-dpng','-r500');
    
    % Plot a rough hitogram
    kk=kk+1;
    figure(kk);clf;
    % edges is unique for each plot
    % but DO NOT affect the data_hist written to file
    % set range of the histogram manually
    switch nn
        case 1
            edges=[115:0.1:123];
        case 2
            edges=[12:0.05:18];
    end
    % set range of the histogram automatically
    bin_size=0.1;
	edges=[min(data_hist):bin_size:max(data_hist)];
    h=histcounts(data_hist,edges,'Normalization', 'probability');
    pp=plot(edges,[0 h],'-o','linewidth',2);
    ax=gca;
    ax.FontSize=20;
%     str_pic=['hist/PLOT_hist-',name1{nn}];
%     print(str_pic,'-dpng','-r500');

    % Write average
    fprintf(id_outfile,'%s\t%.2f (%.2f)\n',[name2{nn}],mean(data_hist),std(data_hist));

    % Write histogram data
    dlmwrite(['hist/DAT_hist-cls1-',name2{nn},'-mdff.dat'],data_hist);
end
fclose(id_outfile);

%     str_pic=['PLOT_ang-rot-ang',num2str(alpha),'-2monomer'];
%     print(str_pic,'-dpng','-r500');