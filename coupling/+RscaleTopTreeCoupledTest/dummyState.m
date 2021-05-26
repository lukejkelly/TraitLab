function state = dummyState(s, prior)
    % Same as makestate but without data and associated calculations
    % Assume catastrophes are already on tree
    state = tree2state(s);

    state.tree = treeclades(state.tree, prior.clade);
    state = UpdateClades(state, [state.leaves, state.nodes], ...
                         size(prior.clade, 2));

    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
    state.length = TreeLength(state.tree, state.root);
    state.kappa = rand;
end
