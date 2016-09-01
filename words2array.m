function array = words2array(wordset,s,L,missingwords)

%function array = words2array(wordset,s)
%the word sets wordset{i} correspond to the nodes s(i)

global LEAF IN OUT MIX

leaves=find([s.type]==LEAF);
NS=length(leaves);
array = OUT.*ones(NS,L);
for k=1:NS
   leaf=leaves(k);
   array(k,wordset{leaf})=IN;
   if nargin==4
       array(k,missingwords{k})=MIX;
   end
end
