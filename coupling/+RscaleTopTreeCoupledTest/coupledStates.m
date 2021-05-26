function [state_x, state_y] = coupledStates(L, theta)
    global ROOT;
    if nargin == 1
        theta = 1e-2;
    end

    tree_x = ExpTree(L, theta);

    numclades = ceil(rand * L / 2);
    clade = synthclades(tree_x, numclades, 2, rand);
    rootmax = 2 * tree_x([tree_x.type] == ROOT).time;
    prior = RscaleTopTreeCoupledTest.clade2prior(clade, rootmax);

    state_x = RscaleTopTreeCoupledTest.dummyState(tree_x, prior);
    state_y = state_x;
end
