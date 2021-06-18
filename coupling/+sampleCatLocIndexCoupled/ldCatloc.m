function ld = ldCatloc(x, state)
    % Uniform over support of catastrophe locations
    catlocs = [state.tree.catloc];
    if any(ismembertol(x, catlocs))
        ld = -log(state.ncat);
    else
        ld = -inf;
    end
end
