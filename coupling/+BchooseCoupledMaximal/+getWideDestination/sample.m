function x = sample(v)
    ind = ceil(length(v) * rand);
    x = v(ind);
end
