function [n,v]=freeprogeny(s,c,i)
 
c1=s(i).child(1);
c2=s(i).child(2);

if ~any(c1==c)
    [n1,v1]=freeprogeny(s,c,c1);
    m1=n1(1);
else
    n1=[];
    m1=s(c1).time;
    v1=[];
end

if ~any(c2==c)
    [n2,v2]=freeprogeny(s,c,c2);
    m2=n2(1);
else
    n2=[];
    m2=s(c2).time;
    v2=[];
end

n=[max(m1,m2),n1,n2];
v=[i,v1,v2];