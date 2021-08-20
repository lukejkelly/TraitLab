function [i, newage, logq] = Schoose(state)
    % Sample node i and new age t_i' ~ U(t_j, t_k) if i not root, and if it is
    % then t_i' ~ U[(t_i + t_j) / 2, 2 * t_i - t_j]
    global ROOT

    s = state.tree;
    i = state.nodes(ceil(rand * (state.NS - 1)));

    iT = s(i).time;
    k = s(i).parent;
    j1 = s(i).child(1);
    j2 = s(i).child(2);
    jT = max(s(j1).time, s(j2).time);

    if s(i).type == ROOT
      %this way to update root should be more adaptive to likelihood
      tau = iT - jT;
      delta = (-0.5 + rand * 1.5) * tau;
      taup = tau + delta;
      logq = log(tau / taup);
      newage = jT + taup;
    else
      newage = jT + rand * (s(k).time - jT);
      logq = 0;
    end

end
