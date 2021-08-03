function t = compareStates(state_x, state_y)
    % Compare state quantities we're interested in (up to relative tolerance)
    fn = fieldnames(state_x);
    ts = zeros(size(fn));
    tv = ones(size(fn));

    for i = 1:length(fn)
        fi = fn{i};
        switch fi
        case {'leaves', 'nodes', 'root', 'cat', 'ncat'}
            ts(i) = all(state_x.(fi) == state_y.(fi));
        case 'claderoot'
            ts(i) = isequaln(state_x.(fi), state_y.(fi));
        case {'mu', 'beta', 'kappa', 'loglkd', 'logprior', 'length'}
            ts(i) = ismembertol(state_x.(fi), state_y.(fi));
        case 'tree'
            ts(i) = equaltrees(state_x.tree, state_y.tree);
        otherwise
            % ignore field
            tv(i) = 0;
        end
    end
    t = all(ts(tv == 1));
end
