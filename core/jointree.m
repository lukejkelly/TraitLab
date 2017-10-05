function s = jointree(sup,sub,ileaf)
% s = jointree(supertree,subtree,ileaf)
% makes one tree out of two where the root of the
% subtree is attached at the node supertree(ileaf)
%GlobalSwitches
global ROOT ANST LEAF

lsup = length(sup);
lsub = length(sub);
% make sure root of sup is second to last node
suproot = find([sup.type]==ROOT);
if suproot ~= lsup - 1
    sup = switchrootposition(sup,suproot,lsup-1);
end
% make sure root of sub is second to last node
subroot = find([sub.type]==ROOT);
if subroot ~= lsub - 1
    sub = switchrootposition(sub,subroot,lsub-1);
    subroot = lsub - 1;
end

s = sup;
% copy supertree adam node to end of new combined tree 
s(lsup+lsub-2) = sup(lsup);
% adjust root node
s(lsup-1).parent = lsup+lsub-2;
% copy subtree nodes into supertree
for i = 1:length(sub)-2
    p = lsup+i-1;
    s(p) = sub(i);
    % adjust numbering
    if s(p).parent ~= subroot
        s(p).parent = s(p).parent + lsup - 1;
    else
        s(p).parent = ileaf;
    end
    s(p).child = s(p).child + lsup - 1;
end    

% put subtree root in correct place
s(ileaf) = sub(subroot);
s(ileaf).parent = sup(ileaf).parent;
s(ileaf).sibling = sup(ileaf).sibling;
% see whether subtree is actually a tree or just offset leaf
if lsub ~=2
    %its a tree
    s(ileaf).child = s(ileaf).child + lsup -1;
    s(ileaf).type = ANST;
else
    % its an offset leaf
    s(ileaf).type = LEAF;
end

% checktree

ok = checktree(s,(lsub+lsup)/2 -1);
if ~isempty(ok)
    ok
    error('Tree is not legal')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = switchrootposition(s,p,np)
global ROOT
% copy node from np
temp = s(np);
% copy node to np
s(np) = s(p);
% fix up parent and child numbering for node at np in correct postion
for i = 1:2
    if s(np).child(i) ~= np
        s(s(np).child(i)).parent = np;
    else
        temp.parent = np;
    end
end
if s(np).type ~= ROOT
    if s(np).parent ~= np
        s(s(np).parent).child(s(np).sibling) = np;
    else
        temp.child(s(np).sibling) = np;
    end
else
    s(s(np).parent).child = np;
end

% copy node saved node into pos p
s(p) = temp;
for i = 1:length(s(p).child)
    s(s(p).child(i)).parent = p;
end
s(s(p).parent).child(s(p).sibling) = p;

ok = checktree(s,length(s)/2);
if ~isempty(ok)
    ok
    error('Tree is not legal')
end
