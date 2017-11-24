clear;
clc;
clf;

prefix = 'cls1-gs3-water';

cccname = [prefix,'_ccc_local.dat'];
cccdata=load(cccname);

rmsdname = [prefix,'_rmsd.dat'];
rmsddata=load(rmsdname);

time_factor = 1e-15*1000*50*1e9;

figure(1);clf;
[ax,p1,p2] = plotyy(cccdata(:,1)*time_factor,cccdata(:,2),rmsddata(:,1)*time_factor,rmsddata(:,2));
set(p1,'color','r');
set(ax(1),'ycolor','r');
ylabel(ax(1),'CCC','FontSize',20);
set(p1,'LineWidth',2);

set(p2,'color','b');
set(ax(2),'ycolor','b');
ylabel(ax(2),'RMSD(A)','FontSize',20);
set(p2,'LineWidth',2);

xlabel('Time (ns)','FontSize',20);
legend({'CCC','RMSD'},'FontSize',20,'Location','SouthEast');

hstr = ['CCC & RMSD profile of MDFF'];
h = title(hstr,'FontSize',20);

pstr = ['PLOT_cls1-gs3-water-ccc-rmsd'];
print(figure(1),'-dpng','-r500',pstr);