function [nstate,U,TOPOLOGY,logq,OK]=RscaleSubTree(state,del,deldel)

%GlobalSwitches;
global LEAF

nstate=state;

nu=progeny(nstate.tree,nstate.root,LEAF);
vec=cumsum(nu(2,:));
tot=vec(2*nstate.NS-1);
r=rand*tot;
i=nu(1,min(find(r<vec)));

if nstate.tree(i).type==LEAF, keyboard;pause; end

if i~=nstate.root
   nu=progeny(nstate.tree,i,LEAF);
end

Vel=nu(1,:);
k=nstate.tree(i).parent; 

rho=del+rand*deldel;

OK=1;

lf=[];
nd=[];
for j=Vel
   if nstate.tree(j).type==LEAF
      lf=[lf,j];
   else
      nd=[nd,j];
   end
end
t0=min([nstate.tree(lf).time]);

for j=nd
   nstate.tree(j).time=t0+rho*(nstate.tree(j).time-t0);
end
for j=lf
   if nstate.tree(nstate.tree(j).parent).time<nstate.tree(j).time
      OK=0;
      break;
   end
end
if nstate.tree(i).time>nstate.tree(k).time, OK=0; end

% Nmovedcats=0;
% 
% for j=[lf,nd];
%     p=state.tree(j).parent;
%     nstate.tree(j).cat=(state.tree(j).cat-state.tree(j).time)/(state.tree(p).time-state.tree(j).time)*(nstate.tree(p).time-nstate.tree(j).time)+nstate.tree(j).time;
%     Nmovedcats=Nmovedcats+length(nstate.tree(j).cat);
% end

if OK
   U=above(nd,nstate.tree,nstate.root);         
   TOPOLOGY=0;
   logq=(length(Vel)-length(lf)-2)*log(rho);
else
   U=[];
   TOPOLOGY=0;
   logq=0;
end

if OK, nstate.length=TreeLength(nstate.tree,nstate.root); end
