clear;
clc;
clf;

prefix = 'cls1-gs3-water_ccc_local';

figure(1);

refname = ['old/result_ccc_local.dat'];
ref=load(refname);
time_factor = 1e-15*1000*10*1e9;
hold on;
plot(ref(1:231,1)*time_factor*10,ref(1:231,2),'color','k','LineWidth',2);
lInfo{1} = ['water(0.1-0.3)'];

cccname = [prefix,'.dat'];
ccc=load(cccname);
time_factor = 1e-15*1000*50*1e9;
plot(ccc(:,1)*time_factor+23,ccc(:,2),'color','r','LineWidth',2);
lInfo{2} = ['water(0.3)'];


xlabel('Time (ns)','FontSize',20);
ylabel('CCC','FontSize',20);
legend(lInfo,'FontSize',20,'Location','SouthEast');

h = title('Overall CCC Profile','FontSize',20);

print(figure(1),'-dpng','-r500','PLOT_cls1-overall-ccc.png');