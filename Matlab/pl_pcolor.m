
function pl_pcolor(x,y)
    N = hist3([x y],'Nbins',[30 30]);
    N_pcolor = N';
    N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
    xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
    yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
    N_pcolor(N_pcolor==0)=NaN;
    h = pcolor(xl,yl,N_pcolor);
    h.FaceColor='interp';
end
