clear;
clc;

str_filename=['cls2-gs3-water-angles-same-row.dat'];
% str_filename=['cls2-gs3-water-angles-diff-row.dat'];

data=load(str_filename);
data(:,1)=data(:,1)/2;

figure(1);clf;
for ii=1:4
    hold on;
    plot(data(:,1),data(:,ii+1),'-','linewidth',2);

    xlabel('Time (ns)','fontsize',20);
    ylabel('Angle (degree)','fontsize',20);
    title('Angle between dimers on same row','fontsize',20);

    str_print=['cls2-gs3-water-angles-same-row'];
end

print(figure(1),'-dpng','-r500',str_print);