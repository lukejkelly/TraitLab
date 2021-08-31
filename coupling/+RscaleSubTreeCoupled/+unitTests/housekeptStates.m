function [state_x, state_y] = housekeptStates(L, theta)
    state_x = RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, theta));
    state_y = housekeeping(state_x, ...
        RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, theta)));
end
