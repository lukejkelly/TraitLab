function [s,c]=ActiveCoversC(s,c,i,LEAF,wk)

if s(i).type>LEAF
   
   c1=s(i).child(1);
   c2=s(i).child(2);

   if s(c1).mark
      [s,c]=ActiveCoversC(s,c,c1,LEAF,wk);
   end
   if s(c2).mark
      [s,c]=ActiveCoversC(s,c,c2,LEAF,wk);
   end
   
   wa=wk;
   wb=wk;
   wc=wk;
   wa(s(c1).clade)=1;
   wb(s(c2).clade)=1;   
   c(wa&wb)=i;
   wd=(wa|wb);
   wc((wa&(~wb))|(wb&(~wa)))=1;
   wc([s(c1).unclade,s(c2).unclade])=1;
   s(i).clade=find(wd);
   s(i).unclade=find(wc);
   
   %s(i).clade=union(s(c1).clade,s(c2).clade);
   %s(i).unclade=unique([setxor(s(c1).clade,s(c2).clade),s(c1).unclade,s(c2).unclade]);
   
   %at the end of this, s(i).clade is the list of clades with at least one
   %member below i
   %s(i).unclade is the list ofclades with at least one member and at least one
   %non-member below i
   %This is input to CleanDownC() where these defs change.
   
else

    c(s(i).clade)=i;
    
end