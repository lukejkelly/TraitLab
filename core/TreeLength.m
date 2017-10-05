function length=TreeLength(s,k)

length=0;
for kc=s(k).child
   length=length+s(k).time-s(kc).time+TreeLength(s,kc);
end
