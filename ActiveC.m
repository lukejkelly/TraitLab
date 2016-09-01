function s=ActiveC(s,i,LEAF,wk)

if s(i).type>LEAF
   
   c1=s(i).child(1);
   c2=s(i).child(2);

   if s(c1).mark==1
      s=ActiveC(s,c1,LEAF,wk);
   end
   if s(c2).mark==1
      s=ActiveC(s,c2,LEAF,wk);
   end
   
   %s(i).clade=[s(c1).clade,s(c2).clade];
   %s(i).unclade=[setxor(s(c1).clade,s(c2).clade),s(c1).unclade,s(c2).unclade];
   
   wa=wk;
   wb=wk;
   wc=wk;
   wa(s(c1).clade)=1;
   wb(s(c2).clade)=1;   
   wd=(wa|wb);
   wc((wa&(~wb))|(wb&(~wa)))=1;
   wc([s(c1).unclade,s(c2).unclade])=1;
   s(i).clade=find(wd);
   s(i).unclade=find(wc);
   
   %s(i).clade=union(s(c1).clade,s(c2).clade);
   %s(i).unclade=unique([setxor(s(c1).clade,s(c2).clade),s(c1).unclade,s(c2).unclade]);
   
end
%