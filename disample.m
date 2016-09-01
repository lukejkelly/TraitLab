function n=disample(p)

%sample n with probability p(n)
c=cumsum(p);
n=1;
r=rand;
while r>c(n)
        n=n+1;
end
