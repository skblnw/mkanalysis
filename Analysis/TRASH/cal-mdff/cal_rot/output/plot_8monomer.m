clear;
clc;

str_filename={'rot_res345-361pa3-barpa1.dat','rot_res345-361pa3-barpa2.dat','rot_res345-361pa3-barpa3.dat'};
str_title={'res345-361pa3-barpa1','res345-361pa3-barpa2','res345-361pa3-barpa3'};

for ii=1:3
    data=load(str_filename{ii});
    data(:,1)=data(:,1)/2;
    switch ii
        case 2
%         data(find(data>90))=180-data(find(data>90));
            data(:,[3 4 7 8])=180-data(:,[3 4 7 8]);
    end
    par_eq=35;
    par_eq=par_eq*20;

    figure(ii);clf;
    cc=hsv(4);
    kk=0;
    data_mean1=mean(data(1,[2:2:end]));
    data_mean2=mean(data(1,[3:2:end]));
    for nn=2:size(data,2)
        kk=int16(nn/2-0.5);
        hold on;
        if rem(nn,2)==0
%             pp(kk)=plot(data(:,1),data(:,nn),'color',cc(kk,:),'linewidth',2);
            plot(data(:,1),data(:,nn),'color',cc(kk,:),'linewidth',2);
            str_legend{nn-1}=['dimer',num2str(kk),'(underneath)'];
            
            data_mean1=[data_mean1 mean(data(par_eq:end,nn))];
        else
            plot(data(:,1),data(:,nn),'--','color',cc(kk,:),'linewidth',2);
            str_legend{nn-1}=['dimer',num2str(kk),'(top)'];
            
            data_mean2=[data_mean2 mean(data(par_eq:end,nn))];
        end
    end
    if ii==1
        legend(str_legend,'fontsize',14,'Location','southeastoutside');
    end
    xlabel('Time (ns)','FontSize',20);
    ylabel('Angle (degree)','FontSize',20);
    xlim([0 50]);
%     ylim([65 95]);
    str_tmp=['Angle (',str_title{ii},') profile'];
    title(str_tmp,'FontSize',16);
    
    str_print=['PLOT_cls2-4dimer-full-rot-',str_title{ii}];
    print(figure(ii),'-r500','-dpng',str_print);
    
    str_outfile=['cls2-4dimer-full-rot-',str_title{ii},'.txt'];
    id_outfile = fopen(str_outfile,'w');
    fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data_mean1(1),mean(data_mean1(2:end)),std(data_mean1(2:end)));
    fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data_mean2(1),mean(data_mean2(2:end)),std(data_mean2(2:end)));
    fclose(id_outfile);
end