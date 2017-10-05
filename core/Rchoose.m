function [s,OK]=Rscale(state,rho)

% GlobalSwitches
[NS,N,L,mu,s,Root,LEAVES,NODES,llkd,olp]=new2old(state);

zt=min([s.time]);

for j=NODES
   s(j).time=zt+(s(j).time-zt)*rho;
end

OK=1;

for k=LEAVES
   if s(s(k).parent).time<s(k).time
      OK=0;
      break;
   end
end

if OK
   state=old2new(NS,L,mu,s,Root,LEAVES,NODES,llkd,olp);
   TOPOLOGY=0;
   [s OK]=MarkRcurs(state,[NODES],TOPOLOGY);
end