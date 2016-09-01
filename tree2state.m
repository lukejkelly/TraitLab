function state=tree2state(s)

global LEAF ROOT ANST

state=pop('state');

state.tree=s;

state.leaves=find([s.type]==LEAF );
state.NS=length([state.leaves]);
state.root=find( [s.type]==ROOT );
state.nodes=[find( [s.type]==ANST ),state.root];
