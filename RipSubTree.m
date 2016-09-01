function [f,r,n]=RipSubTree(s,a)

if s(a).mark
    f=[];
    r=0;
    n=0;
    return;
end

%if s(a).type==LEAF
%    f=TreeNode([],1,[],s(a).time,num2str(1),LEAF);
%    r=1;
%    n=1;
%    return;
%end
    
nc=length(s(a).child);

if nc==0
    f=TreeNode([],1,[],s(a).time,num2str(1),LEAF);
    r=1;
    n=1;
    return;
elseif nc==1
    if s(a).type~=ADAM
        disp('huh? ADAM problems in RipSubTree()');
        keyboard;
    end
    [g.tree,r,m]=RipSubTree(s,s(a).child);
    
    
for k=1:nc
    b=s(a).child(k);
    [g(k).tree,r(k),m(k)]=RipSubTree(s,b,f,n);
end

if nc==1
    if s(n).type~=ADAM
        disp('huh? ADAM problems in RipSubTree()');
    else
    if isempty(g.tree)
        if n~=1
            disp('RipSubTree() built an empty tree - was that unexpected?');
        end
        f(1)=TreeNode([],[],[],realmax,num2str(1),ADAM);
    else
        if n~=(m+1)
            disp('RipSubTree() did not find the expected number of nodes in the subtree');
        else
            f=g.tree;
            f(n)=TreeNode([],[],r,realmax,num2str(n),ADAM);    
            f(r)=TreeNode(n,1,g(r).child,g(r).time,num2str(r),ROOT);
            return;
        end
    end
    

f=[g(1).tree,g(2).tree)];

    nc=length(s(p).child);
for k=1:nc
    cs=s(p).child(k);
    [f,cf,m]=RipSubTree(s,cs,f,n);
    if ~s(p).mark       
        f(n).child(k)=cf;
        f(cf).parent=n;
    end
end

if ~s(p).mark
    f(n).time=s(p).time;  
    cf=n;
    n=n+1;  
end

