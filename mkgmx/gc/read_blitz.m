function A=read_blitz(file);
 %function A=read_blitz(filename);
fid = fopen(file,'r');

%function A=read_blitz(fid)
% A=read_blitz(fid)

str='x'; dim=[];
while (str=='x')
   d=fscanf(fid,'%d',1);    
   str=fscanf(fid,'%s',1);
   dim=[dim d];% disp(str); dim
end;  
assert(str,'[');
if length(dim)<=2 
    A=fscanf(fid,'%g',dim(end:-1:1));
else
    A=zeros(dim([3,2,1]));  %war vorher 2,3,1
    assert(length(dim),3); %geht noch nicht besser
    for i=1:dim(1)
        A(:,:,i)=fscanf(fid,'%g',dim([3,2]));
    end
    %A=read_rec(fid,dim);
end
assert(fscanf(fid,'%s',1),']');

function assert(i,need);
if (i~=need) 
    error(['Format Error in File< character wanted here was',need]);
end
