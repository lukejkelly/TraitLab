function ol = getOverlap(x, y, mcmc)
    % Overlap between distributions of x' ~ p = x * U[a, b] and y' ~ q likewise
    a = mcmc.update.del;
    b = a + mcmc.update.deldel;

    a_p = a * x;
    b_p = b * x;

    a_q = a * y;
    b_q = b * y;

    max_a = max(a_p, a_q);
    min_b = min(b_p, b_q);

    if min_b < max_a
        ol = 0;
    else
        ol = (min_b - max_a) / ((b - a) * max(x, y));
    end
end
