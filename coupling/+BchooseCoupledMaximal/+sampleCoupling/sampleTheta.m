function THETA = sampleTheta()
    % log10(THETA) ~ Uniform(-3, -1)
    THETA = 10^(-3 + rand * 2);
end
