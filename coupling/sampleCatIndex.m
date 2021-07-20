function j = sampleCatIndex(state, i)
    % Sample catastrophe index j uniformly from catastrophes on i
    s = 1:state.cat(i);
    j = discreteUniformSample(s);
end
