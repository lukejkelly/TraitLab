function state = dummyState(s)
    % Same as makestate but without data and associated calculations
    % Assume catastrophes are already on tree
    global LEAF;
    % Randomly offset some leaves if currently isochronous
    if all([s([s.type] == LEAF).time] < 1e-2) && rand < 0.5
        for j = find([s.type] == LEAF)
            s(j).time = rand * (s(s(j).parent).time - s(j).time);
        end
    end
    state = tree2state(s);
    state.claderoot = [];
    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
    state.length = TreeLength(state.tree, state.root);
    state.kappa = rand;
end
