function lp = discreteUniformLogProb(x, s)
    % Log-probability of sampling x uniformly from the entries of s
    n = length(s);
    if n == 0
        error('Support s must be non-empty');
    elseif ismember(x, s)
        lp = -log(n);
    else
        lp = -Inf;
    end
end
