clear;
clc;

data=read_blitz('correl_complex.dat');
figure(1);clf;
plot_corr_matrix(data);
data=read_blitz('correl_rna.dat');
figure(2);clf;
plot_corr_matrix(data);
data=read_blitz('correl_apo.dat');
figure(3);clf;
plot_corr_matrix(data);

for ii=1:3
    figure(ii);
    xtickangle(45);
    pbaspect([1 1 1]);
    ax=gca();
    ax.FontSize=24;
    fig = gcf;
    fig.PaperPositionMode = 'auto'; 
    fig.PaperPosition=[0 0 12 12];
    str_fig=['PLOT_',num2str(ii)];
    print(str_fig,'-dpng');
end