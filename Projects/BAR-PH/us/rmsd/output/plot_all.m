clear;
clc;
clf;

cc=jet(2);
j=0;
%for i=680:15:980
for i=[680 980]
    rmsdname = [ 'us-z',num2str(i),'.rmsd' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*10*1e9
    
    figure(1);
    hold on;
    j=j+1;
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    lInfo{j} = ['z',num2str(i)];
end

xlabel('Time(ns)','FontSize',20);
ylabel('RMSD(A)','FontSize',20);
title('RMSD profile of all 21 windows','FontSize',20);
legend(lInfo,'FontSize',5,'Location','SouthEast');

pstr = ['us-sym-rmsd-time'];
%print(figure(1),'-dpng','-r500',pstr);