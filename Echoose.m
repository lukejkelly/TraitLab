function [i,j,iP,jP,logq,OK]=Echoose(state,mt,prior)
global NARROW OTHER ROOT WIDE

s=state.tree;
Root=state.root;

cladefail=0; %TODO rewrite exchange to find legal moves and choose one

switch mt
    case NARROW
        OK=0;
        N=state.NS-1;
        c=ceil(N*rand);
        iP=state.nodes(c);
        while (iP==state.root)
            c=ceil(N*rand);
            iP=state.nodes(c);
        end
        i=s(iP).child(1+(rand<0.5));
        jP=s(iP).parent;
        j=s(jP).child(OTHER(s(iP).sibling));
        
        if (s(j).time<s(iP).time)
            if (~prior.isclade || (prior.isclade && isequal(s(iP).clade,s(jP).clade)))
                OK=1;
            end
        end
        
        
    case WIDE
        OK=1;
        N=2*state.NS-1;
        while (1)
            i=ceil(N*rand);
            while ( i==Root || ( s(s(i).parent).type==ROOT && s(s(s(i).parent).child(OTHER(s(i).sibling))).time<s(i).time ) ) % i=root or i=(root oldest child)
                i=ceil(N*rand);
            end
            j=ceil(N*rand);
            iP=s(i).parent;
            jP=s(j).parent;
            
            if i~=j && iP~=jP && i~=jP && j~=iP && s(j).time<s(iP).time && s(i).time<s(jP).time
                if (~prior.isclade ||(prior.isclade && isequal(s(iP).clade,s(jP).clade)))
                     break;
                else
                    if cladefail>state.NS, OK=0; break; end %XXX HARDWIRED STOPPING CONDITION
                    cladefail=cladefail+1;
                end
            end
            
        end
        
end

logq=0;