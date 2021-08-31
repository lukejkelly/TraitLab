function [a, b] = getAgePars(state, i, t0, del, deldel)
    % Parameters of uniform distribution
    old = state.tree(i).time;
    a = t0 + del * (old - t0);
    b = t0 + (del + deldel) * (old - t0);
end
