function [state,U,TOPOLOGY]=Eupdate(state,i,j,iP,jP)

% [state,U,TOPOLOGY]=Eupdate(state,i,j,iP,jP)
% switches nodes i and j; iP and jP are the parents of i and j.

s=state.tree;

s(jP).child(s(j).sibling)=i;
s(i).parent=jP;
swap=s(i).sibling;
s(i).sibling=s(j).sibling;
s(iP).child(swap)=j;
s(j).sibling=swap;
s(j).parent=iP;
 
% s(i).cat=(s(i).cat-s(i).time)*(s(jP).time-s(i).time)/(s(iP).time-s(i).time)+s(i).time; % linear transformation
% s(j).cat=(s(j).cat-s(j).time)*(s(iP).time-s(j).time)/(s(jP).time-s(j).time)+s(j).time;

state.tree=s;

U=above([iP,jP],state.tree,state.root);
U=unique(U);
TOPOLOGY=1;

if any([i, j]==state.root), state.length=TreeLength(state.tree,state.root); end