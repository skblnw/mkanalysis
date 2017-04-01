clear;
clc;

str_prefix=['cls1-6tetramer-gs3-contact-I1-'];
figure(1);clf;
cc=hsv(4);
kk=0;
p=0;
for ii=1:4
    kk=kk+1;
    
    str_filename=[str_prefix,num2str(ii),'.dat'];
    data=load(str_filename);
    data(:,1)=data(:,1)/20;
    data(:,2)=data(:,2)/6;
    
    hold on;
    p(kk)=plot(data(:,1),data(:,2),'-o','MarkerSize',3,'color',cc(kk,:),'linewidth',1);

    par_b=50;
    par_B = 1/par_b*ones(1,par_b);
    data_running = filter(par_B,1,data(:,2));
    hold on; plot(data(par_b:end,1),data_running(par_b:end),'color',cc(kk,:),'LineWidth',4)

    str_legend{kk}=['pair#',num2str(ii)];
end
legend(p,str_legend,'fontsize',12,'location','northeastoutside');
xlabel('Time (ns)','FontSize',20);
ylabel('# of contacts','FontSize',20);
str_tmp=['Contact profile of Region I1'];
title(str_tmp,'FontSize',20);

str_print=['PLOT_cls1-6tetramer-gs3-contact-I1'];
print(figure(1),'-r500','-dpng',str_print);

str_prefix=['cls1-6tetramer-gs3-contact-I2-'];
figure(2);clf;
cc=hsv(6);
kk=0;
p=0;
for ii=1:6
    kk=kk+1;
    
    str_filename=[str_prefix,num2str(ii),'.dat'];
    data=load(str_filename);
    data(:,1)=data(:,1)/20;
    data(:,2)=data(:,2)/6;
    
    hold on;
    p(kk)=plot(data(:,1),data(:,2),'-o','MarkerSize',3,'color',cc(kk,:),'linewidth',1);

    par_b=50;
    par_B = 1/par_b*ones(1,par_b);
    data_running = filter(par_B,1,data(:,2));
    hold on; plot(data(par_b:end,1),data_running(par_b:end),'color',cc(kk,:),'LineWidth',4)

    str_legend{kk}=['pair#',num2str(ii)];
end
legend(p,str_legend,'fontsize',12,'location','northeastoutside');
xlabel('Time (ns)','FontSize',20);
ylabel('# of contacts','FontSize',20);
str_tmp=['Contact profile of Region I2'];
title(str_tmp,'FontSize',20);

str_print=['PLOT_cls1-6tetramer-gs3-contact-I2'];
print(figure(2),'-r500','-dpng',str_print);