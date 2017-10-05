function [state OK]=MarkRcurs(state,VERTICES,TOPOLOGY, ignoreearlywarn)

global za1 za2 zb1 zb2 zq1 zq2;

if nargin<=3, ignoreearlywarn=0; end

OK=1;

if numel(VERTICES)>0
    [state.tree(VERTICES).mark]=deal(1);
end

if TOPOLOGY
   za1=zeros(1,state.L);
   za2=zeros(1,state.L);
   zb1=zeros(1,state.L);
   zb2=zeros(1,state.L);
   zq1=zeros(1,state.L);
   zq2=zeros(1,state.L);
   [state.tree]=ActiveI(state.tree,state.root);
   state.tree(state.root).CovI=[state.tree(state.root).ActI{:}];
   %[state.tree,state.claderoot]=CoversI(state.tree,state.root,state.claderoot);
   [state.tree]=CoversI(state.tree,state.root);
end

[state.tree toolarge]=WUprune(state.tree,state.root,state.mu,state.p,state.kappa,state.cat);

if toolarge % This happens when the value of mu is too large. 
% The llkd is basically -infinity, so we want to reject, but such values
% lead to bugs in WUprune. RJR 15/06/11
    OK=0;
    if ~ignoreearlywarn
        disp('Proposed value of death rate is too large, leading to numerical issues. If this happens a lot, you might have convergence issues. If this only happens at the beginning of the run, you should be fine.');
    end
end

if numel(VERTICES)>0
    [state.tree(VERTICES).mark]=deal(0);
end
