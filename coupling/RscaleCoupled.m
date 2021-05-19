function [var_x, var_y] = RscaleCoupled(state_x, state_y, a, b)
    % Function which attempts to maximally couple root times and returns
    % corresponding variation terms for inputting to Rscale

    i = state_x.root;
    if i ~= state_y.root
        error('Root indices do not match');
    end

    t0_x = min([state_x.tree(state_x.leaves).time]);
    t0_y = min([state_y.tree(state_y.leaves).time]);

    old_x = state_x.tree(i).time;
    old_y = state_y.tree(i).time;

    % new ~ min + U(min + a * (old - min), min + b * (old - min))
    a_p = t0_x + a * (old_x - t0_x);
    b_p = t0_x + b * (old_x - t0_x);

    a_q = t0_y + a * (old_y - t0_y);
    b_q = t0_y + b * (old_y - t0_y);

    [new_x, new_y] = maximalCouplingUniform(a_p, b_p, a_q, b_q);

    var_x = (new_x - t0_x) / (old_x - t0_x);
    var_y = (new_y - t0_y) / (old_y - t0_y);
end
