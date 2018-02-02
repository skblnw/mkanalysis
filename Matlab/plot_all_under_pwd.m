clear;
clc;
close all;

% Import sort_nat
addpath('../../../../../../scripts/Plot/Matlab/plugins');

files=dir('*.dat');
filenames={files.name};
filenames=sort_nat(filenames);

str_tmp=['TABLE_dat.txt'];
id_outfile = fopen(str_tmp,'w');
kk=0;

% Time
% either set ns_start or frame_start
psperframe=1000;
f2ns=psperframe/1000;
frame_start=300;
% ns_start=0;

for jj=filenames
    data_output=[];
    data_bar=[];
    data_err=[];
    
    % Get name of the file
    tmp1=strsplit(char(jj),'.');
    tmp2=strsplit(char(tmp1(1)),'-');
    name=char(tmp2(2));
    data=load(char(jj));
    data(:,1)=data(:,1)*f2ns;

    % If plot in separate figures
    if false
        kk=kk+1;
        figure(kk);clf;
    else
        kk=kk+1;
        figure(1);
    end
    
    for ii=2:length(data(1,:))
        % If plot running averages
        if true
            cc=hsv(length(files));
            hold on;
            p_dat(kk)=plot(data(:,1),data(:,2),'-','color',cc(kk,:),'LineWidth',1);
            p_dat(kk).Color(4)=0.2; % Transparency

            par_b=10;
            par_B = 1/par_b*ones(1,par_b);
            data_running = filter(par_B,1,data(:,2));
            hold on;
            p_run(kk)=plot(data(par_b:end,1),data_running(par_b:end),'color',cc(kk,:),'LineWidth',3);
            
            str_legend{kk}=['resid ',name];
        else
            hold on;
            plot(data(:,1),data(:,ii),'--','linewidth',1);
        end
    end

    % If multiple columns
    if false
        for ii=1:length(data)
            data2(ii,1)=data(ii,1);
            data2(ii,2)=mean(data(ii,2:end));
        end
        p_avg=plot(data2(:,1),data2(:,2),'-','linewidth',2);
    end


    legend(p_run,str_legend,'fontsize',16, 'location','bestoutside');
%     title(name, 'fontsize',20);
%     xlabel('ns', 'fontsize',20);
%     ylabel(str_y{jj}, 'fontsize',20);

    ax=gca;
    ax.FontSize=20;
%     xlim([0 1000]);
%     ax.XTick=[0:50:250];
%     ylim([0 10]);
%     ax.YTick=[0:2:7];

%     print(['PLOT_',name],'-dpng','-r100');
    fprintf(id_outfile,'%s\t%.2f\t(%.2f)\n',name,mean(data(frame_start:end,2)),std(data(frame_start:end,2)));
end
fclose(id_outfile);