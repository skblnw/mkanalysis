clear;
clc;

par_threshold=0;
str_filename={'contact-i1.dat','contact-i2-under.dat','contact-i2-above.dat'};
str_title={'Region I1','Region I2 (above)','Region I2 (underneath)'};
for nn=1:3
    data=load(str_filename{nn});
    data(:,1)=data(:,1)/2;

    figure(nn);clf;
    
    kk=0;
    for ii=1:361
        if max(data(:,ii+1))>par_threshold
            kk=kk+1;
        end
    end
    
    cc=hsv(kk);
    kk=0;
    p=0;
    for ii=1:361
        if max(data(:,ii+1))>par_threshold
            kk=kk+1;
            hold on;
            p(kk)=plot(data(:,1),data(:,ii+1),'--o','MarkerSize',3,'color',cc(kk,:),'linewidth',.5);

            par_b=50;
            par_B = 1/par_b*ones(1,par_b);
            data_running = filter(par_B,1,data(:,ii+1));
            hold on; plot(data(par_b:end,1),data_running(par_b:end),'color',cc(kk,:),'LineWidth',3)

            str_legend{kk}=['resid',num2str(ii)];
        end
    end
%     h_leg=legend(p,str_legend,'location','eastoutside');
    gridLegend(p,3,str_legend,'fontsize',10,'location','eastoutside');
    xlabel('Time (ns)','FontSize',20);
    ylabel('# of contacts','FontSize',20);
    str_tmp=['Contact profile of ',str_title{nn},''];
    title(str_tmp,'FontSize',20);
    
    str_print=['PLOT_cls2-4dimer-full-contact-',num2str(nn),''];
    print(figure(nn),'-r500','-dpng',str_print);
end