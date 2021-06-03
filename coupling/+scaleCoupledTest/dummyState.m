function state = dummyState()
    state = tree2state(ExpTree(5, 5e-2));
    state.lambda = rand;
    state.mu = rand;
    state.beta = rand;
    state.kappa = 0.25 + rand * 0.75;
    state.rho = rand;
    state.nu = state.lambda * state.kappa / state.mu;
    state.length = TreeLength(state.tree, state.root);
    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
end
