function [state_x, state_y] = coupledStates(L, theta)
    [state_x, state_y] = deal(...
        RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, theta)));
end
