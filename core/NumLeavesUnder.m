function m=NumLeavesUnder(a,s)

global LEAF

n=progeny(s,a,LEAF);
m=sum(n(2,:)==0);