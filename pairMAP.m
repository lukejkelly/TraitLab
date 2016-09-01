function mode=pairMAP(a,b,mu,p)

A=length(intersect(a,b));
L=length(union(a,b));
if A>0;
   b=A/(L+A);
   c=sqrt(3*p^2-2*p+1+2*L*p^2/A);
   x=(b/p^2)*(p-1+c);
   mode=-log(x)/mu;
else
   mode=inf;
end