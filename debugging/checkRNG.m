function checkRNG(n)
    % store RNG state, print n draws from rand, restore the RNG seed
    if nargin == 0
        n = 1;
    end
    s = rng;
    fprintf([repmat('%g ', 1, n - 1), '%g\n'], rand(1, n));
    rng(s);
end
