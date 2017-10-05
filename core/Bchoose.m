function [i,j,k,newage,logq]=Bchoose(state,mt,THETA,prior)

global ROOT WIDE OTHER
N=2*state.NS-1;
s=state.tree;

FAIL=0;

i=ceil(N*rand);
while s(i).type==ROOT
   i=ceil(N*rand);
end
iP=s(i).parent;
iT=s(i).time;

if mt==WIDE
    
    if prior.isclade, 
        r=[]; 
        for a=1:2*state.NS, 
            if s(a).type<ROOT & length(s(a).clade)>=length(s(iP).clade) & length(s(s(a).parent).clade)>=length(s(iP).clade) & isempty(setdiff(s(iP).clade,s(a).clade)) & isempty(setdiff(s(s(a).parent).clade,s(iP).clade)), 
                r=[r,a]; 
            end; 
        end
        N=length(r); 
    else
        r=1:N;
    end
    
    if N>4
        j=r(ceil(N*rand));
        k=s(j).parent;
        while ( s(k).time<=iT || i==j || i==k )
            j=r(ceil(N*rand));
            k=s(j).parent;
        end
    else
        FAIL=1;
        k=-1;
        j=-1; 
    end
    
else
   
   if s(iP).type==ROOT
      FAIL=1;
      k=-1;
      j=-1;
   else
      k=s(iP).parent;
      j=s(k).child(OTHER(s(iP).sibling));
   end
   
end

if FAIL || k==iP || j==iP
   
   newage=[];
   logq=-inf;
   
elseif s(j).type==ROOT
   
   jT=s(j).time;
   delta=-(1/THETA)*log(rand);                            %technical, check
   newage=jT+delta;
   
   PiP=s(iP).parent;
   CiP=s(iP).child(OTHER(s(i).sibling));
   
   PiPT=s(PiP).time;
   CiPT=s(CiP).time;
   old_minage=max(iT,CiPT);
   old_range=PiPT-old_minage;
   
   q=exp(delta*THETA)/(THETA*old_range);                 %XXX technical, check
   logq=log(q);
   
elseif s(iP).type==ROOT
   
   jT=s(j).time;
   kT=s(k).time;
   new_minage=max(iT,jT);
   new_range=kT-new_minage;
   newage=new_minage+rand*new_range;
   
   CiP=s(iP).child(OTHER(s(i).sibling));
   CiPT=s(CiP).time;
   q=exp((CiPT-s(iP).time)*THETA)*new_range*THETA;     %XXX technical, check
   logq=log(q);
      
else
   
   jT=s(j).time;
   kT=s(k).time;
   new_minage=max(iT,jT);
   new_range=kT-new_minage;
   newage=new_minage+rand*new_range;
   
   PiP=s(iP).parent;
   CiP=s(iP).child(OTHER(s(i).sibling));
   
   PiPT=s(PiP).time;
   CiPT=s(CiP).time;
   old_minage=max(iT,CiPT);
   old_range=PiPT-old_minage;
   
   q=new_range/old_range;
   logq=log(q);
   
end
