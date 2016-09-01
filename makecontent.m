function content = words2content(true)

GlobalSwitches;

content = pop('content');

content.array = OUT.*ones(true.NS,true.L);
for k=1:true.state.NS
   leaf=true.state.leaves(k);
   content.array(k,true.wordset{leaf})=IN;
end

content.NS=true.state.NS;
content.L=true.state.L;
content.cognate=true.cognate;
content.language=true.language;
