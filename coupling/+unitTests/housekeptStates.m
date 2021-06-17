function [state_x, state_y] = housekeptStates(L, theta)
    state_x = unitTests.dummyState(ExpTree(L, theta));
    state_y = housekeeping(state_x, unitTests.dummyState(ExpTree(L, theta)));
end
