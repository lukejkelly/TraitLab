function s=swap(s,j,k)

if j==k
    return;
end

%order so only special case is j has k as parent
if ~isempty(s(k).parent) & s(k).parent==j
   x=j;j=k;k=x;
end

so=s;

sj=s(j);
sk=s(k);

s(j)=sk;
s(k)=sj;

if isempty(sj.parent) | sj.parent~=k
    %'ne'
    if ~isempty(sj.parent)     
        s(sj.parent).child(sj.sibling)=k;
    end
    if ~isempty(sj.child) 
        [s(sj.child).parent]=deal(k);
    end
    if ~isempty(sk.parent)
        s(sk.parent).child(sk.sibling)=j;
    end
    if ~isempty(sk.child)  
        [s(sk.child).parent]=deal(j);
    end
    %result=checktree(s,size(so,2)/2);
    %if ~isempty(result), disp('error pp'); keyboard;pause; end
else
    %'e'
    if ~isempty(sj.child)
        [s(sj.child).parent]=deal(k);
    end
    s(k).parent=j;
    s(j).parent=sk.parent;
    if ~isempty(sk.parent)
        s(sk.parent).child(sk.sibling)=j;
        s(j).sibling=sk.sibling;
    end
    c=sk.child(sk.child~=j);
    s(j).child=[k,c];
    s(k).sibling=1;
    if ~isempty(c), s(c).sibling=2; s(c).parent=j; end
    %result=checktree(s,size(so,2)/2);
    %if ~isempty(result), disp('error pp'); keyboard;pause; end
end