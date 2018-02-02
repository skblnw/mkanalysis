% /*
% ################################
% ## MATLAB script to plot 3D vectors from a matrix of 6 col using quiver3
% ##   Most original version of quiver3 MATLAB code series, mostly has been
% ##   replaced by quiver1 (1 frame only), quiverc (cumulative), quiverm
% ##   (multiple models), etc.
% ## Kevin Apr 2014
% ## Input: vec.dat (6 col matrix)
% ## Output: Simple MATLAB plot or even a GIF if you uncomment lines after
% ##         "drawnow"
% ## Units:
% ## Notes: the input file must include coordinates which is inconvenient
% ################################
% */

filename='vec.dat';
outname='sym-fn-100-hold.gif';
figure
for i=1:40:4000
    k=i+39;
    dlmstr=['A',num2str(i),'..F',num2str(k)];
    data=dlmread(filename,' ',dlmstr);
    %figure
    quiver3(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),'b');
    fn=k/40;
    tstr=['sym-Frame ',num2str(fn)];
    title(tstr);
    hold on
    axis ([10 70 -20 20 -10 70])
    view(-2,8)
    hold off
    
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
% 	[A,map] = rgb2ind(im,256); 
%         if i == 1;
%             imwrite(A,map,outname,'gif','LoopCount',Inf,'DelayTime',0.1);
%         else
%             imwrite(A,map,outname,'gif','WriteMode','append','DelayTime',0.1);
%         end
end
