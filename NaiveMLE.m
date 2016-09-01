function nMLE=NaiveMLE(a,b,mu)

li=length(intersect(a,b));
if li>0;
   x=2*li/(length(a)+length(b));
   nMLE=-0.5*log(x)/mu;
else
   nMLE=inf;
end