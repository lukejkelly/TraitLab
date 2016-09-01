function v=catperlength(s,k,v)

for kc=s(k).child
   v(kc,1:2)=[length(s(kc).cat),(s(k).time-s(kc).time)];
   v=catperlength(s,kc,v);
end
