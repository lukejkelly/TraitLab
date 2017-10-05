function [new, q1, q2]=GetLegal(s,old,root)
% [new, q1, q2]=GetLegal(state,old)
% Given tree s and node old, returns a neighbour of old that is allowed to
% hold a catastrophe. Also returns q1 and q2, the number of available neighbours
% for old and new.
% RJR, 07/05/07

global OTHER

poss=setdiff([s(old).child s(old).parent s(s(old).parent).child(OTHER(s(old).sibling))],root);
q1=length(poss);
if q1==0, keyboard; disp('Error in Getlegal.m: node does not have a sibling nor children and its parent is the root.'); end
new=poss(ceil(rand*q1));
q2=length(setdiff([s(new).child s(new).parent s(s(new).parent).child(OTHER(s(new).sibling))],root));