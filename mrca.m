function p=mrca(U,s,Root)

% function p=mrca(U,s,Root)
%
% Find the MRCA of a set U of leaf nodes on tree s
%
% Seems to work for sets U including ancestral node labels unless 
% any node in the set U is a descendent of any other node or U=[Root]
% CHANGE GKN 18/10/02

global ROOT

if nargin==2
   Root=find([s.type]==ROOT);
end

for p=U
   s(p).mark=1;
   while p~=Root
      p=s(p).parent;
      if s(p).mark==1
         break;
      else
         s(p).mark=1;
      end
   end
end

p=Root;
FINISHED=0;
while ~FINISHED
   c=[s(p).child];
   if isempty(c) | sum([s(c).mark])==2
      FINISHED=1;
   else
      p=c([s(c).mark]==1);
   end
end