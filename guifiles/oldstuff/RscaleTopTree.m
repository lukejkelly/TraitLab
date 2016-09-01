function [nstate,U,TOPOLOGY,logq,OK]=RscaleTopTree(state,prior,del,deldel)

%GlobalSwitches;


nstate=state;

[n,v]=freeprogeny(state.tree,union(state.leaves,state.claderoot(prior.upboundclade)),state.root);
 
rho=del+rand*deldel;

OK=1;

check=v;
for j=1:length(v)
   d=min([state.tree(state.leaves).time]);
   nstate.tree(v(j)).time=d+rho*(state.tree(v(j)).time-d);
   check=[check,state.tree(v(j)).child];
end
check=unique(check);
for j=check
   if nstate.tree(state.tree(j).parent).time<nstate.tree(j).time
      OK=0;
      break;
   end
end
 
% Nmovedcats=0;
% 
% for j=check;
%     p=state.tree(j).parent;
%     nstate.tree(j).cat=(state.tree(j).cat-state.tree(j).time)/(state.tree(p).time-state.tree(j).time)*(nstate.tree(p).time-nstate.tree(j).time)+nstate.tree(j).time;
%     Nmovedcats=Nmovedcats+length(nstate.tree(j).cat;
% end


if OK
   U=above(check,nstate.tree,nstate.root);         
   TOPOLOGY=0;
   logq=(length(v)-2)*log(rho);
else
    %keyboard;
    %disp('TopTree proposed a bad tree');
   U=[];
   TOPOLOGY=0;
   logq=0;
end

if OK, nstate.length=TreeLength(nstate.tree,nstate.root); end
