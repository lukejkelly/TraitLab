function s=CleanDownC(s,i,LEAF,wk)

if s(i).type>LEAF
    
    c1=s(i).child(1);
    c2=s(i).child(2);
    
    if s(c1).mark
        s=CleanDownC(s,c1,LEAF,wk);
    end
    
    if s(c2).mark
        s=CleanDownC(s,c2,LEAF,wk);
    end
    
    wa=wk;
    wb=~wk;
    wa(s(i).clade)=1;
    wb(s(i).unclade)=0;
    s(i).clade=find(wa&wb);
    
end