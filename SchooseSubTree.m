function [i,k,Vel]=SchooseSubTree(state)

GlobalSwitches;

nu=progeny(state.tree,state.root,LEAF);
vec=cumsum(nu(2,:));
tot=vec(2*state.NS-1);
r=rand*tot;
i=nu(1,min(find(r<vec)));

if state.tree(i).type==LEAF, keyboard;pause; end

if i~=state.root
   nu=progeny(state.tree,i,LEAF);
end

Vel=nu(1,:);
k=state.tree(i).parent; 
