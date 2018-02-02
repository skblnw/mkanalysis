clear;
clc;
prefix={'output_z/'};
sel=[2];
name={'zaxis'};

str_tmp=['TABLE_avg.txt'];
id_outfile = fopen(str_tmp,'w');
str_tab=[''];
kk=0;

% Time
% either set ns_start or frame_start
psperframe=1000;
f2ns=psperframe/1000;
ns_eq=10;
frame_eq=ns_eq*1000/psperframe;

for jj=1:1
    for nn=1:sel(jj)
        filename = [prefix{jj},'sel',num2str(nn),'.dat'];
        data=load(filename);
        data(:,1)=data(:,1)*f2ns;
        data_avg=[];
        data_sum=[];

		% Plot multiple columns
        kk=kk+1;
        figure(kk);clf;
        for ii=2:size(data,2)
            hold on;
            plot((1:size(data,1))*f2ns,data(:,ii),'-','linewidth',2);
        end

        for ii=1:size(data,1)
            data_avg(ii)=mean(data(ii,2:end));
            data_sum(ii)=sum(data(ii,2:end));
        end
%         p=plot((1:size(data_avg,2))*f2ns,data_avg,'-','linewidth',2);
%         plot((1:size(data_sum,2))*f2ns,data_sum,'-','linewidth',2);

		% Write mean
        fprintf(id_outfile,'%s\t%.2f (%.2f)\n',str_tab,mean(data_avg(end-frame_eq:end)),std(data_avg(end-frame_eq:end)));

        legend(name,'fontsize',20, 'location','best');
        title('', 'fontsize',20);
        xlabel('ns', 'fontsize',20);
%         ylabel('', 'fontsize',20);

        ax=gca;
        ax.FontSize=20;
        if nn==1
%             ylim([20 120]);
        end
%         xlim([]);
%         ax.XTick=[0:50:250];
%         ax.YTick=[0:2:7];

        str_pic=['PLOT_',name{jj},'-sel',num2str(nn)];
%         print(str_pic,'-dpng');
		
    end
end
fclose(id_outfile);