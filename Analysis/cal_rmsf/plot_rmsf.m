clear;
clc;

name={'PROT','MONO','BAR','PH'};
folder={'run','run1','run2'};
ps=[50 100 100];
leg={'sim1','sim3','sim2'};

figure(1);clf;
for jj=[3]
    for ii=4:4
        for kk=1:2
            filename = [folder{jj},'/rmsf_P',num2str(kk),'-',name{ii},'.dat'];
            data=load(filename);
            hold on;
            plot(data(:,1),data(:,2), 'linewidth',2);
        end
    end
end
 
legend({'protein'},'fontsize',20, 'location','best');
% title('Ups1', 'fontsize',20);
xlabel('Residue', 'fontsize',20);
ylabel('RMSF(A)', 'fontsize',20);

ax=gca;
ax.FontSize=20;
% ax.XTick=[0:20:180];
% ax.YTick=[0:2:7];
 
% print('PLOT_rmsf','-dpng','-r500');