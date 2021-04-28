function [var_x, var_y] = RscaleRoot(state_x, state_y, a, b)
    % Function which attempts to maximally couple root times and returns
    % corresponding variation terms, where new_time = var * old_time, for
    % inputting to Rscale

    i = state_x.root;
    if i ~= state_y.root
        error('Root indices do not match');
    end

    old_x = state_x.tree(i).time;
    old_y = state_y.tree(i).time;

    % new ~ U(old * a, old * b)
    [new_x, new_y] = maximalCouplingUniformScaling(old_x, old_y, a, b);

    var_x = new_x / old_x;
    var_y = new_y / old_y;
end
