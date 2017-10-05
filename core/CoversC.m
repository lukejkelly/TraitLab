function [s,c]=CoversC(s,i,c,LEAF,wk)

c1=s(i).child(1);
c2=s(i).child(2);

if s(c1).mark 
   if s(c1).type>LEAF
      [s,c]=CoversC(s,c1,c,LEAF,wk);
   else
      c(s(c1).clade)=c1;
   end
end

if s(c2).mark 
   if s(c2).type>LEAF
      [s,c]=CoversC(s,c2,c,LEAF,wk);
   else
      c(s(c2).clade)=c2;
   end
end

%c([s(c1).clade,s(c2).clade])=i;
wa=wk;
wb=wk;
wa(s(c1).clade)=1;
wb(s(c2).clade)=1;
c(wa&wb)=i;
