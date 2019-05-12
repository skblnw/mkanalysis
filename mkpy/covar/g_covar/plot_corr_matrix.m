function h=plot_corr_matrix(A,bsmooth);
%plot_corr_matrix(A,[bsmooth])
%plots a matrix in a colored panel. 
%if bsmooth>0 it will smooth the results in both dimensions
%to have slow color changes.
if (nargin<2)
        bsmooth=0;
end
R=1:size(A,1);
D=A;
if (bsmooth)
B=expand_matrix(D(R,R),5);
H=gaussav2d(B,6,1.5);
L=5;
else
H=D; 
L=1;
end
nH=zeros(size(H,1)+1,size(H,2)+1);
nH(1:size(H,1),1:size(H,2))=H;
h=pcolor((1:size(nH,1))/L-.5,(1:size(nH,2))/L-.5,nH);%,'EdgeColor','none');
set(h,'EdgeColor','None');
if (bsmooth)
  set(h,'FaceColor','Interp');
end
axis([0.5 length(R)-0.5 0.5 length(R)-0.5]);
step=ceil(length(R)/20);
set(gca,'Xtick',1:step:length(R));
set(gca,'XtickLabel',num2str(R(1:step:end)'),'TickDir','out');
set(gca,'Ytick',1:step:length(R));
set(gca,'YtickLabel',num2str(R(1:step:end)'),'TickDir','out');
end

function [H]=gaussav2d(data,av_points,sigma)
%hist2d_av2(M,data,av_points,sigma)
%macht 2d histogramm der daten in den ersten 
% beiden Spalten von data
% wenn av_points <>0 dann
% wird jeder datenpunkt nach einem Gaussgewichtung
% auf mehrere umliegende Gitterpunkte verteilt
% M is number of gridpoints in each dimension
     

    H=zeros(size(data,1)+av_points*2,size(data,2)+av_points*2);
    Gh=mache_2D_gewicht(av_points,sigma);
    Gh=Gh./sum(sum(Gh));
    for x=1:size(data,1)
        for y=1:size(data,2)
           for k=-av_points:av_points
               for l=-av_points:av_points
                   H(x+k+av_points,y+l+av_points)=H(x+k+av_points,y+l+av_points)...
                       +Gh(k+1+av_points,l+1+av_points)*data(x,y);
               end
           end
        end
    end     
    H=H(1+av_points:end-av_points,1+av_points:end-av_points);       
end

function Gh=mache_2D_gewicht(M,sigma);
   [X,Y]=meshgrid(-M:M,-M:M);
   Gh=1./(2*pi)/sigma^2.*exp(-X.^2/2/sigma^2-Y.^2/2/sigma^2);
end

function A=expand_matrix(B,E);
%expand_matrix(matrix,expand)
%expands matrix by coping every data point to ExE elements
%if you smooth with gaussav2d later it might even look good
A=zeros(size(B,1)*E,size(B,2)*E);
for i=1:size(B,1)
    for j=1:size(B,2)
        for k=1:E
            for l=1:E
                A((i-1)*E+k,(j-1)*E+l)=B(i,j);
            end
        end
    end
end
end

function [H]=gaussav2(data,av_points,sigma)
%hist2d_av2(M,data,av_points,sigma)
%macht 2d histogramm der daten in den ersten 
% beiden Spalten von data
% wenn av_points <>0 dann
% wird jeder datenpunkt nach einem Gaussgewichtung
% auf mehrere umliegende Gitterpunkte verteilt
% M is number of gridpoints in each dimension
     

    H=zeros(size(data,1)+av_points*2,size(data,2)+av_points*2);
    Gh=mache_2D_gewicht(av_points,sigma);
    Gh=Gh./sum(sum(Gh));
    for x=1:size(data,1)
        for y=1:size(data,2)
           for k=-av_points:av_points
               for l=-av_points:av_points
                   H(x+k+av_points,y+l+av_points)=H(x+k+av_points,y+l+av_points)...
                       +Gh(k+1+av_points,l+1+av_points)*data(x,y);
               end
           end
        end
    end     
    H=H(1+av_points:end-av_points,1+av_points:end-av_points);       
end
      