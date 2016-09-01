function n=progeny(s,i,LEAF)

%Returns a list of all nodes below or equal to i, and the number of descendants each node has (excluding itself).

if s(i).type==LEAF %isempty(s(i).child)
  n=[i;0];
else
  c1=s(i).child(1);
  c2=s(i).child(2);
  n1=progeny(s,c1,LEAF);
  n2=progeny(s,c2,LEAF);
  n=[ [i;2+n1(2,1)+n2(2,1)], n1, n2 ];
end
