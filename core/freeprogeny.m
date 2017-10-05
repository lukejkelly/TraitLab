function [n,v]=freeprogeny(s,c,d,i,prior)
 
%at the moment this assumes (line 11) a leaf must be fixed, or at least in an
%upper-bounded range. Could easilly have leaves with completely unknown
%ages - this should be corrected in the calling program by exclusing those
%leaves from c - but also needs work here.

c1=s(i).child(1);
c2=s(i).child(2);

if ~any(c1==c)
    [n1,v1]=freeprogeny(s,c,d,c1,prior);
    m1=n1(1);
else
    j=find(c1==d);
    if isempty(j)
        m1=s(c1).time;
    else
         if length(j)==1
            m1=prior.clade{prior.upboundclade(j)}.lowlim;
        else
            wka=[prior.clade{prior.upboundclade(j)}];
            m1=max([wka.lowlim]);
        end
    end
    n1=[];
    v1=[];
end

if ~any(c2==c)
    [n2,v2]=freeprogeny(s,c,d,c2,prior);
    m2=n2(1);
else
    j=find(c2==d);
    if isempty(j)
        m2=s(c2).time;
    else
        if length(j)==1
            m2=prior.clade{prior.upboundclade(j)}.lowlim;
        else
            wka=[prior.clade{prior.upboundclade(j)}];
            m2=max([wka.lowlim]);
        end
    end
    n2=[];
    v2=[];
end

n=[max(m1,m2),n1,n2];
v=[i,v1,v2];