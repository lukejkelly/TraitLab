function t = checkCoupling(state_x, state_y)
    % If state variables match then attempt to compare trees
    t = checkCoupling.compareStates(state_x, state_y) ...
        && checkCoupling.compareTrees(state_x, state_y);

    % TODO: change after further checks
    if t && ~checkCoupling.extraChecks(state_x, state_y)
        warning('Some state fields do not match');
    end
end
