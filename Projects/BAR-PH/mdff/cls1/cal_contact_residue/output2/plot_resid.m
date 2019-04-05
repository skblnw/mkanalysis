clear;
clc;

str_filename={'contact-i1.dat','contact-i2-above.dat','contact-i2-under.dat'};
cons_start=120;
data_plot=[];
data_std=[];
for nn=1:3
    data_raw=load(str_filename{nn});
    data_raw=data_raw/2;

    for ii=2:size(data_raw,2)
        jj=ii-1;
        if mean(data_raw(cons_start:end,ii)) < 0.0
            data_tmp(jj)=0;
            data_tmpsd(jj)=0;
        else
            data_tmp(jj)=mean(data_raw(cons_start:end,ii));
            data_tmpsd(jj)=std(data_raw(cons_start:end,ii));
        end
    end

    data_plot=[data_plot;data_tmp];
    data_std=[data_std;data_tmpsd];
end

for nn=1:3
    data_index=[];
    data_find=[];
    data_findsd=[];
    for ii=find(data_plot(nn,:))
        data_index=[data_index ii];
        data_find=[data_find data_plot(nn,ii)];
        data_findsd=[data_findsd data_std(nn,ii)];
    end
    data_index=[data_index;data_find;data_findsd];

    [tmp_Y,tmp_I]=sort(data_index(2,:),'descend');
    data_find=data_index(:,tmp_I);
    str_outfilename=['data_find',num2str(nn),'.txt'];
    id_outfile = fopen(str_outfilename,'w');
    fprintf(id_outfile,'%3.0f\t%3.1f(%0.1f)\n',data_find);
    fclose(id_outfile);
end

figure(1);clf;
bb = bar(data_plot');
% bb = barwitherr(data_std',1:361,data_plot');
set(bb(1),'facecolor','g');
set(bb(1),'edgecolor','g');
set(bb(2),'facecolor','r');
set(bb(2),'edgecolor','r');
set(bb(3),'facecolor','b');
set(bb(3),'edgecolor','b');

xlim([0 361]);
set(gca,'XTick',1:20:361);
xlabel('Residue ID','FontSize',20);
ylabel('# of contacts','FontSize',20);
legend({'Region I1','Region I2 (above)','Region I2 (underneath)'},'FontSize',16,'Location','northwest');
% legend('boxoff');
title('# of contacts at residues','FontSize',24);
grid on;

str_print=['PLOT_cls1-4dimer-full-resseq-all'];
% print(figure(1),'-r500','-dpng',str_print);


figure(2);clf;
par_range={[1 17],[18 97],[98 111],[112 133],[134 148],[149 179],[180 233],[234 249],[250 270],[271 290],[291 361]};
for nn=1:length(par_range)
    subplot(length(par_range)/1,1,nn);
    bb = barwitherr(data_std',1:361,data_plot');
    
    set(bb(1),'facecolor','g');
    set(bb(1),'edgecolor','g');
    set(bb(2),'facecolor','r');
    set(bb(2),'edgecolor','r');
    set(bb(3),'facecolor','b');
    set(bb(3),'edgecolor','b');

    xlim(par_range{nn});
    if ismember(nn,[1 3 4 6 8 9 10])
        ylim([0 15]);
        set(gca,'YTick',[5 10]);
        set(gca,'fontsize',6);
        set(gca,'XTick',par_range{nn}(1):1:par_range{nn}(2));
        grid on;
    else
        set(gca,'YTick',[]);
        set(gca,'fontsize',6);
        set(gca,'XTick',par_range{nn}(1):10:par_range{nn}(2));
    end
end

set(gcf,'NextPlot','add');
axes;
h = title('# of contacts at residues in details','FontSize',24);
set(gca,'Visible','off');
set(h,'Visible','on');

str_print=['PLOT_cls1-4dimer-full-resseq-detail'];
print(figure(2),'-r500','-dpng',str_print);