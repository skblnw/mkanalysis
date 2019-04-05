clc;
clear;

prefix1 = ['old/symmetry-2-1-32'];
prefix2 = ['cls1-gs3-water'];

% Plot bar charts
bdata1=[];
bdata2=[];
bdata3=[];
for i=1:4
    string=[prefix1,'-contact-I1-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/4;
    
    bdata1 = [bdata1 idata(1,2)];
    
    string=[prefix2,'-contact-I1-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/4;
    bdata3 = [bdata3 mean(idata(:,2))];
end

for i=2:5
    string=[prefix1,'-contact-I2-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/4;
    
    bdata1 = [bdata1 idata(1,2)];
    
    string=[prefix2,'-contact-I2-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/4;
    bdata3 = [bdata3 mean(idata(:,2))];
end

figure(1);clf;
b = bar([bdata1;bdata3]');
set(b(1),'facecolor','k');
set(b(2),'facecolor','r');

set(gca, 'XTick',1:2, 'XTickLabel',{'' ''});
ylabel('# of contacts','FontSize',20);
legend({'initial','water(0.3)'},'FontSize',20,'Location','northwest');
title('Statistics comparing initial and final structures after MDFF','FontSize',18);
print(figure(1),'-dpng','-r500','PLOT_cls1-gs3-water_contact_stat.png');