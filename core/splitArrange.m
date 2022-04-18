function x=splitArrange(s,n,a,b,i)
% vertical coordinates of nodes for draw.m
% renamed from split.m to avoid collision with built-in function

global LEAF

if s(i).type==LEAF %isempty(s(i).child)
  x=[i;a+(b-a)/2];
else
  c1=s(i).child(1);
  c2=s(i).child(2);
  if (n(c1)+n(c2))==0
    x=[[i;(a+b)/2],splitArrange(s,n,a,(a+b)/2,c1),splitArrange(s,n,(a+b)/2,b,c2)];
  else
    dnt=(b-a)/(n(c1)+n(c2));
    if n(c1)==0
     %d1=dnt/2;
     d1=dnt;
    elseif n(c2)==0
       %d1=(n(c1)-1/2)*dnt;
     d1=(n(c1)-1)*dnt;
    else
     d1=n(c1)*dnt;
    end
    if n(c1)>n(c2)
      x=[[i;(b-d1+(a+b)/2)/2],splitArrange(s,n,b-d1,b,c1),splitArrange(s,n,a,b-d1,c2)];
    else
      x=[[i;(a+d1+(a+b)/2)/2],splitArrange(s,n,a,a+d1,c1),splitArrange(s,n,a+d1,b,c2)];
    end
  end
end
