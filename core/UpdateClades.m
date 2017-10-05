function state=UpdateClades(state,VERTICES,nc)

%GlobalSwitches;
global LEAF

[state.tree(VERTICES).mark]=deal(1);
work=zeros(1,nc);
%[state.tree]=ActiveC(state.tree,state.root,LEAF,work);
%[state.tree,state.claderoot]=CoversC(state.tree,state.root,state.claderoot,LEAF,work);
[state.tree,state.claderoot]=ActiveCoversC(state.tree,state.claderoot,state.root,LEAF,work);
[state.tree]=CleanDownC(state.tree,state.root,LEAF,work);

[state.tree(VERTICES).mark]=deal(0);
