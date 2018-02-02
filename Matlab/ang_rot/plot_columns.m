close all;
clear;
clc;
name={'output/ang-rot-PH-BAR-'};
str_y={'Angle (^{\circ})'};

kk=0;
str_tmp=['TABLE_ang-rot.txt'];
id_outfile = fopen(str_tmp,'w');
for jj=1:3
    data_output=[];
    data_bar=[];
    data_err=[];
    for nn=1:2
        filename = [name{1},num2str(jj),'.dat'];
        data=load(filename);
        f2ns=200;
        data(:,1)=data(:,1)*f2ns/1000;
        
        if nn==1
            data=[data(:,1) data(:,2:2:end)];
        else
            data=[data(:,1) data(:,3:2:end)];
        end

        switch jj
        case 2
            data_tmp=data;
            data_tmp(find(data_tmp>90))=180-data_tmp(find(data_tmp>90));
            data=data_tmp;
        end
        
        kk=kk+1;
        figure(kk);clf;
        for ii=2:length(data(1,:))
            hold on;
            plot(data(:,1),data(:,ii),'--','linewidth',1);
        end

        for ii=1:length(data)
            data_avg(ii,1)=data(ii,1);
            data_avg(ii,2)=mean(data(ii,2:end));
        end
        p=plot(data_avg(:,1),data_avg(:,2),'-','linewidth',2);

        legend(p,'avg.','fontsize',20, 'location','best');
        title('', 'fontsize',20);
        xlabel('ns', 'fontsize',20);
        ylabel(str_y{1}, 'fontsize',20);

        ax=gca;
        ax.FontSize=20;
        % ax.XTick=[0:50:250];
        % ax.YTick=[0:2:7];

        str_pic=['PLOT_ang-rot-ang',num2str(jj),'-',num2str(nn)];
        print(str_pic,'-dpng','-r100');
        str_tab=['ang',num2str(jj),'-',num2str(nn)];
        fprintf(id_outfile,'%s\t%.2f (%.2f)\n',str_tab,mean(data_avg(1000:end,2)),std(data_avg(1000:end,2)));
        
        data3(:,1)=data_avg(:,1);
        data3(:,nn+1)=data_avg(:,2);
    end
    kk=kk+1;
    figure(kk);clf;
    hold on;
    plot(data3(:,1),data3(:,2),'b-','linewidth',2);
    plot(data3(:,1),data3(:,3),'r-','linewidth',2);
    
    switch jj
    case 1
        legend({'P-M','P-P'},'fontsize',40, 'location','northeast','orientation','horizontal');
    end
%     title('', 'fontsize',20);
%     xlabel('ns', 'fontsize',20);
%     ylabel(str_y{1}, 'fontsize',20);

    ax=gca;
    ax.FontSize=64;
%     xlim([0 30]);
    switch jj
    case 1
        ylim([0 40]);
        ax.YTick=[0 20 40];
        ax.XTick=[];
    case 2
        ylim([50 110]);
        ax.YTick=[50 80 110];
        ax.XTick=[];
    case 3
        ylim([70 120]);
        ax.YTick=[70 95 120];
        ax.XTick=[0 15 30];
    end
    str_pic=['PLOT_ang-rot-ang',num2str(jj),'-2monomer'];
%     print(str_pic,'-dpng','-r500');
end
fclose(id_outfile);