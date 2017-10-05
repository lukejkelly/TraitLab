function [NS,N,L,mu,s,Root,LEAVES,NODES,claderoot,llkd,olp]=new2old(state)

%function [NS,N,L,mu,s,Root,LEAVES,NODES,claderoot,llkd,olp]=new2old(state);
%
%Temporary (?) : deal the contents of the new 'state' data structure
%into the old variable names

s=state.tree;
mu=state.mu;
Root=state.root;
llkd=state.loglkd;
LEAVES=state.leaves;
NODES=state.nodes;
L=state.L;
NS=state.NS;
N=2*NS-1;
olp=state.logprior;
claderoot=state.claderoot;