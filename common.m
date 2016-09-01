function p=common(i,j,s,Root)

% function p=common(i,j,s,Root)
%
% Use mrca([i,j],s,Root)
% mrca() is more general, and only very slightly slower than common()
% CHANGE GKN 18/10/02

p=i;
s(p).mark=1;
while p~=Root
   p=s(p).parent;
   s(p).mark=1;
end

p=j;
while s(p).mark==0
   p=s(p).parent;
end