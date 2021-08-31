function ld = catlocLogProb(x, state, i)
    % Uniform over support of catastrophe locations on branch i
    if any(ismembertol(x, state.tree(i).catloc))
        ld = -log(state.cat(i));
    else
        ld = -inf;
    end
end
