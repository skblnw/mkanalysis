clc;
clear;

% prefix1 = ['old/symmetry-2-1-23'];
prefix2 = ['cls2-gs3-water'];
par_eq=50;
par_eq=par_eq*20;

% Plot bar charts
bdata1=[];
bdata2=[];
bdata2std=[];
for i=1:4
    string=[prefix2,'-contact-I1-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/2;
    
    bdata1 = [bdata1 idata(1,2)];
    bdata2 = [bdata2 mean(idata(par_eq:end,2))];
    bdata2std = [bdata2std std(idata(par_eq:end,2))];
end

for i=1:6
    string=[prefix2,'-contact-I2-',num2str(i),'.dat'];
    idata=importdata(string);
    idata(:,2)=idata(:,2)/2;
    
    bdata1 = [bdata1 idata(1,2)];
    bdata2 = [bdata2 mean(idata(par_eq:end,2))];
    bdata2std = [bdata2std std(idata(par_eq:end,2))];
end

figure(1);clf;
b = bar([bdata1;bdata2]');
set(b(1),'facecolor','k');
set(b(2),'facecolor','r');

set(gca, 'XTick',1:2, 'XTickLabel',{'' ''});
ylim([0 90]);
ylabel('# of contacts','FontSize',20);
legend({'initial','water(0.3)'},'FontSize',20,'Location','northwest');
title('Statistics comparing contact residue pairs before and after MDFF','FontSize',16);
print(figure(1),'-dpng','-r500','PLOT_cls2-gs3-water_contact_stat.png');

str_outfile=['DAT_cls2-gs3-water_contact_stat.txt'];
id_outfile = fopen(str_outfile,'w');
fprintf(id_outfile,'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',bdata1);
for ii=1:10
    fprintf(id_outfile,'%.2f (%.2f)\t',bdata2(ii),bdata2std(ii));
end
fclose(id_outfile);