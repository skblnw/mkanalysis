clear;
clc;

str_prefix=['output2/contact-A'];

figure(1);clf;
cc=hsv(11);
kk=0;
res=[];

str_tmp=['TABLE_contacts-anyother','.txt'];
id_outfile = fopen(str_tmp,'w');

frame_start=60;

for nn=3:13
    kk=kk+1;
    str_name=[str_prefix,num2str(nn),'-anyother.dat'];
    data=load(str_name);

    pp(kk)=plot(data(:,1),data(:,2),'-','color',cc(kk,:),'LineWidth',1);
	pp(kk).Color(4)=0.2; % Transparency

%     par_b=10;
%     par_B = 1/par_b*ones(1,par_b);
%     data_running = filter(par_B,1,data(:,2));
%     hold on;
%     plot(data(par_b:end,1),data_running(par_b:end),'color',cc(kk,:),'LineWidth',3);

    data_smooth = moving_average(data(:,2),100);
    hold on;
    pp(kk)=plot(data(:,1),data_smooth,'color',cc(kk,:),'LineWidth',3);

    str_legend{kk}=['A',num2str(nn)];

    fprintf(id_outfile,'A%d\t%.2f\t%.2f(%.2f)\n',nn,data(1,2),mean(data(frame_start:end,2)),std(data(frame_start:end,2)));
    res=[res mean(data(frame_start:end,2))];
    
    grid on;
    xlabel('Time (ns)','fontsize',20);
    ylabel('Contact #','fontsize',20);
%     str_tmp=['Contacts ',str_title{nn}];
%     title(str_tmp,'fontsize',20);
    legend(pp,str_legend,'fontsize',20,'location','bestoutside');
    ax=gca;
    ax.FontSize=20;
    ax.XTick=[0:50:250];
end

str_tmp=['PLOT_contacts-anyother'];
% print(figure(1),'-depsc','-loose',str_tmp);

fprintf(id_outfile,'Avg.\t%.2f(%.2f)\n',mean(res),std(res));
fclose(id_outfile);