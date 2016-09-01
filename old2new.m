function state=old2new(NS,L,mu,s,Root,LEAVES,NODES,claderoot,OLL,olp)

state=pop('state');

state.NS=NS;
state.L=L;
state.mu=mu;
state.tree=s;
state.root=Root;
state.leaves=LEAVES;
state.nodes=NODES;
state.loglkd=OLL;
state.logprior=olp;
state.claderoot=claderoot;
