function x = discreteUniformSample(s)
    % Sample x uniformly from the entries of s
    n = length(s);
    i = ceil(rand * n);
    x = s(i);
end
