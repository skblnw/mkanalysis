clear;
clc;

str_filename1={'cls1-gs3-water-angles-diff-row.dat','cls1-gs3-water-angles-same-row.dat'};
str_filename2={'cls2-gs3-water-angles-diff-row.dat','cls2-gs3-water-angles-same-row.dat'};
str_title={'ang-diff-row','ang-same-row'};

data1=load(str_filename1{1});
data2=load(str_filename2{1});
data1(:,1)=data1(:,1)/2;
data2(:,1)=data2(:,1)/2;
par_eq1=25;
par_eq1=par_eq1*20;
par_eq2=35;
par_eq2=par_eq2*20;


figure(1);clf;
plot(data1(:,1),data1(:,2),'-r','linewidth',2);
xlabel('Time (ns)','fontsize',20);
ylabel('Angle (degree)','fontsize',20);
% hold on;
figure(2);clf;
plot(data2(:,1),data2(:,2),'-g','linewidth',2);
xlabel('Time (ns)','fontsize',20);
ylabel('Angle (degree)','fontsize',20);

% legend({'cls1','cls2'},'fontsize',20,'location','east');


str_tmp=['Dimers in different row']; 
% title(str_tmp,'FontSize',24);

str_print=['PLOT_ang-bar-bar'];
print(figure(1),'-dpng','-r500','PLOT_ang-bar-bar-mdff1');
print(figure(2),'-dpng','-r500','PLOT_ang-bar-bar-mdff2');
% print(figure(1),'-painters','-depsc','-loose',str_print);

str_outfile=['compare-full-',str_title{1},'.txt'];
id_outfile = fopen(str_outfile,'w');
fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data1(1,2),mean(data1(par_eq1:end,2)),std(data1(par_eq1:end,2)));
fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data2(1,2),mean(data2(par_eq2:end,2)),std(data2(par_eq2:end,2)));
fclose(id_outfile);

data1=load(str_filename1{2});
data2=load(str_filename2{2});
data1(:,1)=data1(:,1)/2;
data2(:,1)=data2(:,1)/2;
data1(:,2)=(data1(:,2)+data1(:,3))/2;
data2(:,2)=(data2(:,2)+data2(:,3))/2;
data1(:,3)=(data1(:,4)+data1(:,5))/2;
data2(:,3)=(data2(:,4)+data2(:,5))/2;
for ii=1:1
    figure(3);clf;
    plot(data1(:,1),data1(:,ii+1),'-r','linewidth',2);
    xlabel('Time (ns)','fontsize',20);
    ylabel('Angle (degree)','fontsize',20);
    leg_str{1}=['cls1'];
%     hold on;
    figure(4);clf;
    plot(data2(:,1),data2(:,ii+1),'-g','linewidth',2);
    leg_str{2}=['cls2'];
end
% legend(leg_str,'fontsize',20,'location','southeast');
xlabel('Time (ns)','fontsize',20);
ylabel('Angle (degree)','fontsize',20);
str_tmp=['Dimers in the same row'];
% title(str_tmp,'FontSize',24);
str_print=['PLOT_compare-full-',str_title{2}];
print(figure(3),'-dpng','-r500','PLOT_ang-ph-bar-mdff1');
print(figure(4),'-dpng','-r500','PLOT_ang-ph-bar-mdff2');
% print(figure(2),'-painters','-depsc','-loose',str_print);

str_outfile=['compare-full-',str_title{2},'.txt'];
id_outfile = fopen(str_outfile,'w');
fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data1(1,2),mean(data1(par_eq1:end,2)),std(data1(par_eq1:end,2)));
fprintf(id_outfile,'%.2f\t%.2f(%.2f)\n',data2(1,2),mean(data2(par_eq2:end,2)),std(data2(par_eq2:end,2)));
fclose(id_outfile);