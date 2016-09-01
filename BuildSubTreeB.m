function f=BuildSubTreeB(s,mask,language)

global ADAM

[s.mark]=deal(0);
VERTICES=[mask,[state.tree(mask).parent]];
[state.tree(VERTICES).mark]=deal(1);

a=find(s.type==ADAM);
f=RipSubTree(s,a);