function L=MyPoisson(alpha);

L=0;
tot=0;
while tot<=alpha,
   tot=tot-log(rand);
   L=L+1;
end
L=L-1;