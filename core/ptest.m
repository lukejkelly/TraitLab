function Y=ptest(x0,n)

Y=zeros(1,n);
x=x0;

for k=1:(n-1)
   Y(k)=x;
   a=0.5+1.5*rand;
   xp=x^a;
   MHR=xp/(a*x);
   if rand<min(1,MHR)
      x=xp;
   end
end
Y(n)=x;