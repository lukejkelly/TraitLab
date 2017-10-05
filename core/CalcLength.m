function l=CalcLength(state)

global ANST ROOT

l=0;
for node=1:length(state.tree)
    if ((state.tree(node).type == ANST) | (state.tree(node).type == ROOT))
        c1=state.tree(node).child(1);
        c2=state.tree(node).child(2);
        l=l + 2*state.tree(node).time - state.tree(c1).time - state.tree(c2).time;
    end
end