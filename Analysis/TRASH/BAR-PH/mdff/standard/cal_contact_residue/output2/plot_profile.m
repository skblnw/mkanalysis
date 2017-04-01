clear;
clc;

par_scale=2;
par_threshold=2;
str_filename={'contact-i1.dat','contact-i2-above.dat','contact-i2-under.dat'};
str_title={'Region I1','Region I2 (above)','Region I2 (underneath)'};
for nn=1:3
    data=load(str_filename{nn});
    data(:,1)=data(:,1)/2;

    figure(nn);clf;
    
    kk=0;
    for i=2:size(data,2)
        if max(data(:,i))>par_threshold
            kk=kk+1;
        end
    end
    
    cc=hsv(kk);
    kk=0;
    p=0;
    for ii=2:size(data,2)
        if max(data(:,ii))>par_threshold
            kk=kk+1;
            hold on;
            p(kk)=plot(data(:,1),data(:,ii)/par_scale,'-o','MarkerSize',3,'color',cc(kk,:),'linewidth',1);

            par_b=50;
            par_B = 1/par_b*ones(1,par_b);
            data_running = filter(par_B,1,data(:,ii));
            hold on; plot(data(par_b:end,1),data_running(par_b:end)/par_scale,'color',cc(kk,:),'LineWidth',4)

            str_legend{kk}=['resid',num2str(ii-1)];
        end
    end
    legend(p,str_legend,'fontsize',12,'location','northeastoutside');
    xlabel('Time (ns)','FontSize',20);
    ylabel('# of contacts','FontSize',20);
    str_tmp=['Contact profile of ',str_title{nn}];
    title(str_tmp,'FontSize',20);
    
    str_print=['PLOT_cls1-4dimer-full-contact-',num2str(nn)];
    print(figure(nn),'-r500','-dpng',str_print);
end