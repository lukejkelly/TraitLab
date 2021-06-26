function d = pg_cdf(x, p)
    if x >= 0
        d = catsOnBranches.pg_pdf(x, p) + catsOnBranches.pg_cdf(x - 1, p);
    else
        d = 0;
    end
end
