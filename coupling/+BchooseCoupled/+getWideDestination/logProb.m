function l = logProb(x, v)
    if ismember(x, v)
        l = -log(length(v));
    else
        l = -Inf;
    end
end
