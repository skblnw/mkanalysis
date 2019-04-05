clear;
clc;
clf;

figure(1);clf;
cc=hsv(4);
j=1;
% for i=[616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980]
for i=[616 630 640 650]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(2);clf;
cc=hsv(4);
j=1;
for i=[660 670 680 695]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(3);clf;
cc=hsv(4);
j=1;
for i=[710 725 740 755]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(4);clf;
cc=hsv(4);
j=1;
for i=[765 770 785 800]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(5);clf;
cc=hsv(4);
j=1;
for i=[815 830 845 860]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(6);clf;
cc=hsv(4);
j=1;
for i=[875 890 905 920]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

figure(7);clf;
cc=hsv(4);
j=1;
for i=[935 950 965 980]
    rmsdname = [ 'us-z',num2str(i),'_rmsd.dat' ];
    rmsd=load(rmsdname);
    
    time_factor=2e-15*1000*100*1e9;
    
    subplot(2,2,j);
    plot(rmsd(:,1)*time_factor,rmsd(:,2),'color',cc(j,:),'LineWidth',2);
    
    xlabel('Time(ns)','FontSize',20);
    ylabel('RMSD(A)','FontSize',20);
%     title('RMSD profile of all 28 windows','FontSize',20);
    str_leg=['z',num2str(i)];
    h_leg=legend(str_leg);
    set(h_leg,'FontSize',20,'Location','SouthEast');

    j=j+1;
end

for i=1:7
    pstr=['us-sym-rmsd-',num2str(i)];
    print(figure(i),'-dpng','-r500',pstr);
end