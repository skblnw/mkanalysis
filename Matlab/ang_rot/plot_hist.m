close all;
clear;
clc;
name={'output1/ang-rot-PH-BAR-'};
str_y={'Angle (^{\circ})'};

kk=0;
str_tmp=['hist/TABLE_ang-rot.txt'];
id_outfile = fopen(str_tmp,'w');
for alpha=1:3
    data_output=[];
    data_bar=[];
    data_err=[];
    for nn=1:2
        filename = [name{1},num2str(alpha),'.dat'];
        data=load(filename);
        f2ns=10;
        data(:,1)=data(:,1)*f2ns/1000;
        
        if nn==1
            data=[data(1000:end,2:2:end)];
        else
            data=[data(1000:end,3:2:end)];
        end

        switch alpha
            case 1
                edges=[0:40];
            case 2
                data_tmp=data;
                data_tmp(find(data_tmp>90))=180-data_tmp(find(data_tmp>90));
                data=data_tmp;
                edges=[50:110];
            case 3
                edges=[70:120];
            case 4
                data_tmp=data(:,2:end);
                data_tmp(find(data_tmp>90))=180-data_tmp(find(data_tmp>90));
                data(:,2:end)=data_tmp;
        end
        
        data2=reshape(data,[],1);
        
        kk=kk+1;
        figure(kk);clf;
        h=histcounts(data2,edges,'Normalization', 'probability');
        plot(edges,[0 h]);
        
        str_tab=['alpha',num2str(alpha),'-angle',num2str(nn)];
        fprintf(id_outfile,'%s\t%.2f (%.2f)\n',str_tab,mean(data2),std(data2));

        dlmwrite(['hist/DAT_cls1-alpha',num2str(alpha),'-angle',num2str(nn),'.dat'],data2);
    end
%     str_pic=['PLOT_ang-rot-ang',num2str(alpha),'-2monomer'];
%     print(str_pic,'-dpng','-r500');
end
fclose(id_outfile);