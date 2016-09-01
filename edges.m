function r=edges(s,i,leafnames)

%edges compute splits of tree s, return as splits x leaves binary matrix

els=s(s(i).parent).time-s(i).time;
if isempty(s(i).child)
   r=zeros(1,length(leafnames));
   i=find(strcmp(s(i).Name,leafnames));
   r(i)=1;
   r=[els,r];
else
   c1=s(i).child(1);
   c2=s(i).child(2);
   r1=edges(s,c1,leafnames);
   r2=edges(s,c2,leafnames);
   r=[[els,(sum([r1(:,2:end);r2(:,2:end)])>eps)];r1;r2];
end