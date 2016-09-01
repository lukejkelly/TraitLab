function B=equaltrees(s,t)

global ANST ROOT;

N=length(s);
M=length(t);

if N==M
    B=1;
    r=find([t.type]==ROOT);
    for i=1:N
        if s(i).type==ANST
            vn={s(listleaves(s,i)).Name};
            v=GetLeaves(t,vn);
            if ~IsClade(v,t,r), B=0; return; end
        end
    end
    
end

end

function v=listleaves(s,i)

global LEAF;

if s(i).type==LEAF
    v=i;
else
    v=[listleaves(s,s(i).child(1)),listleaves(s,s(i).child(2))];
end

end