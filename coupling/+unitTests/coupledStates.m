function [state_x, state_y] = coupledStates(L, theta)
    [state_x, state_y] = deal(unitTests.dummyState(ExpTree(L, theta)));
end
