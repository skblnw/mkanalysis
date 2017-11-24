clear;
clc;
clf;

prefix = 'cls2-gs3-water';

cccname = [prefix,'_ccc_local.dat'];
cccdata=load(cccname);

rmsdname = [prefix,'_rmsd.dat'];
rmsddata=load(rmsdname);

time_factor = 1e-15*1000*50*1e9;

figure(1);clf;
plot(cccdata(:,1)*time_factor,cccdata(:,2),'color','r','LineWidth',2);
ylabel('CCC','FontSize',20);

xlabel('Time (ns)','FontSize',20);
legend({'CCC'},'FontSize',20,'Location','SouthEast');

hstr = ['CCC profile of MDFF'];
h = title(hstr,'FontSize',20);

pstr = ['PLOT_cls1-gs3-water-ccc'];
print(figure(1),'-dpng','-r500',pstr);

figure(2);clf;
plot(rmsddata(:,1)*time_factor,rmsddata(:,2),'color','b','LineWidth',2);
ylabel('CCC','FontSize',20);

xlabel('Time (ns)','FontSize',20);
legend({'RMSD'},'FontSize',20,'Location','SouthEast');

hstr = ['RMSD profile of MDFF'];
h = title(hstr,'FontSize',20);

pstr = ['PLOT_cls1-gs3-water-rmsd'];
print(figure(2),'-dpng','-r500',pstr);
