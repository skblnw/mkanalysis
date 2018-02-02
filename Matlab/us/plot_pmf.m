clear;
clc;

figure(1);clf;
cc=cool(5);
kk=0;

min_range=7;
max_range=7.5;

for ii=[40]
    filename = ['output_W',num2str(ii),'_M0/bsResult-',num2str(ii),'ns-1000.xvg'];
    fileID = fopen(filename);
    formatSpec = '%f %f %f';
    data = textscan(fileID, formatSpec, 'CommentStyle', {'#','@'}, 'HeaderLines', 16);
    fclose(fileID);
    
    min_pmf=min(data{2});
    min_index=find(data{2}==min_pmf);
    min_sep=data{1}(min_index);
    str_tmp=['TABLE_pmfmin.txt'];
    id_outfile = fopen(str_tmp,'w');
    fprintf(id_outfile,'%s\t%d\t%f\n',filename,min_index,min_sep);
    fclose(id_outfile);
    
    ref_index=find(data{1}>min_range & data{1}<max_range);
    ref_pmf=data{2}(ref_index);
    ref_avg=mean(ref_pmf);
    data{2}=data{2}-ref_avg;
    
    kk=kk+1;
    hold on;
%     tmp=shadedErrorBar(data{1},data{2},data{3},{'linewidth',4,'color','black'},1);
%     tmp=struct2cell(tmp);
%     pp(kk)=tmp{1};
    plot(data{1},data{2},'linewidth',4,'color','black');
    str_leg{kk}=['Last ',num2str(ii),'ns'];
end
xlim([min(data{1}) max(data{1})]);

if 0
total=70;
window=30;
move=10;
for ii=1:5
    filename = ['output_W',num2str(window),'_M',num2str(move),'/pmf-',num2str(ii),'-1000.xvg'];
    fileID = fopen(filename);
    formatSpec = '%f %f';
    data = textscan(fileID, formatSpec, 'CommentStyle', {'#','@'}, 'HeaderLines', 15);
    
    ref_index=find(data{1}>min_range & data{1}<max_range);
    ref_pmf=data{2}(ref_index);
    ref_avg=mean(ref_pmf);
    data{2}=data{2}-ref_avg;
    
    kk=kk+1;
    hold on;
    pp(kk)=plot(data{1},data{2},'linewidth',1,'color',cc(kk-1,:));
    str_leg{kk}=[num2str(total-move*(ii-1)-window),'-',num2str(total-move*(ii-1)),'ns'];
end
end

legend(str_leg,'fontsize',24,'location','southeast');
title('PMF Profile','fontsize',20);
xlabel('Separation (nm)','fontsize',20);
ylabel('PMF (kcal/mol)','fontsize',20);

grid on;
ax=gca;
ax.XGrid='off';
ax.FontSize=20;
% xlim([min(data{1}) max(data{1})]);
% ylim([0 8]);
% ax.YTick=[0:1:8];

% print(['FIG_pmf_compare_block'],'-dpng');