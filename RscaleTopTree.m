function [nstate,U,TOPOLOGY,logq,OK]=RscaleTopTree(state,prior,del,deldel)

nstate=state;

%GKN change 22/3/11 to add treatment of adamrange clades (to rootrange
%clades)
if ~prior.isclade
    tops=[];
    warning('TraitLab:MCMC','RscaleTopTree update called without clades - should use regular scale.');
else
    tops=state.claderoot(prior.upboundclade);
    if any(prior.isupboundadamclade) %upbound adam clades
        ubac=prior.upboundclade(prior.isupboundadamclade);
        tops(prior.isupboundadamclade)=[state.tree(state.claderoot(ubac)).parent];
    end
end    
[n,v]=freeprogeny(state.tree,[state.leaves,tops],tops,state.root,prior);
     
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
