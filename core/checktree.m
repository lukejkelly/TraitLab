function result=checktree(s,NS)

global ROOT ADAM ANST LEAF

%%
%[x,y]=draw(s,102,1,'newtree');drawnow;      
    N=2*NS;
    result=[];    
    Root=find([s.type]==ROOT);
    if NS~=(size(s,2)/2)
        result=[result,[0.1;-1]]
        return;
    end
%%
    for k=1:N
        if s(k).type==ROOT
            if k~=Root
                result=[result,[1.1;k]]
            end
            if s(s(k).parent).type~=ADAM
                result=[result,[1.2;k]]
            end
            if s(s(k).parent).child~=Root
                result=[result,[1.3;k]]
            end
        elseif s(k).type==LEAF
            if ~isempty(s(k).child)
                result=[result,[2.2;k]]
            end
            if s(s(k).parent).child(s(k).sibling)~=k
                result=[result,[2;k]]
            end
            if s(s(k).parent).time<s(k).time
                result=[result,[2.1;k]]
            end
        elseif s(k).type==ADAM
            %if k~=2*NS
            %    result=[result,[4.4;k]]
            %end
            if ~isempty(s(k).parent)
                result=[result,[4.1;k]]
            end
            if size(s(k).child,2)~=1 || s(k).child~=Root
                result=[result,[4.2;k]]
            end
            if s(k).time<s(Root).time
                result=[result,[4.3;k]]
            end
        else
            if s(k).type~=ANST
                result=[result,[3.1;k]]
            end
            if s(s(k).parent).child(s(k).sibling)~=k
                result=[result,[3.2;k]]
            end
            if s(s(k).parent).time<s(k).time
                result=[result,[3.3;k]]
            end
        end
    end