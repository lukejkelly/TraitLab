function [mObs, mExp, vObs, vExp] = getMoments(d, n)
    [~, ~, k] = LogRhoPrior(0);
    b = 1 / k;

    D = sum(d, 2);
    N = sum(n, 2);
    P = b / (b + D);

    mObs = mean(N);
    mExp = catsOnBranches.pg_mu(P);
    vObs = var(N);
    vExp = catsOnBranches.pg_var(P);
end
