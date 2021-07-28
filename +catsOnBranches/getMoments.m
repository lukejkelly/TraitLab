function [mObs, mExp, vObs, vExp] = getMoments(d, n)
    % Moments assuming samples are independent
    [~, ~, k] = LogRhoPrior(0);
    b = 1 / k;

    D = sum(d, 2);
    N = sum(n, 2);
    P = b ./ (b + D);

    % E[N]
    mObs = mean(N);
    mExpGivenTree = arrayfun(@(p) catsOnBranches.pg_mu(p), P);
    mExp = mean(mExpGivenTree);

    % V[N]
    vObs = var(N);
    vExpGivenTree = arrayfun(@(p) catsOnBranches.pg_var(p), P);
    vExp = mean(vExpGivenTree) + var(mExpGivenTree);
end
