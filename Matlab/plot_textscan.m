clear;
clc;
 
j=1;
cc=hsv();
for i=[]
    filename = ['',num2str(i), ''];
    fileID = fopen(filename);
    formatSpec = '%f %f';
    data = textscan(fileID, formatSpec, 'CommentStyle', {'#' ,'@' }, 'HeaderLines', 15);
   
    figure(1);
    plot(data{1},data{2}, 'linewidth',2,'color' ,cc(j,:));
    xlim([]);
    str_leg{j}=[''];
    hold on;
    j=j+1;
end
 
legend(str_leg,'fontsize',20, 'location','southeast' );
title('', 'fontsize',20);
xlabel('', 'fontsize',20);
ylabel('', 'fontsize',20);
 
pstr = [''];
print(figure(1),'-dpng','-r500',pstr);
