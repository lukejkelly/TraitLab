function b=below(s,i)

% function b=below(s,i)
%
% Make a list of nodes below i on tree (not including i)
% GKN 28/10
% TODO tidy up duplication in this, progeny.m, SchooseSubTree.m and so on

b=[];
for k=s(i).child
  b=[b,k,below(s,k)];
end

