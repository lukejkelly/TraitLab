function t = compareTrees(state_x, state_y)
    % Node types and cat counts already compared at state level
    [s_x, l, n] = deal(state_x.tree, state_x.leaves, state_x.nodes);
    s_y = state_y.tree;

    fn = fieldnames(s_x);
    ts = zeros(size(fn));
    tv = ones(size(fn));

    for i = 1:length(ts)
        fi = fn{i};
        switch fi
        case 'Name'
            ts(i) = all(strcmp({s_x(l).Name}, {s_y(l).Name}));
        case 'xi'
            ts(i) = all(arrayfun(@ismembertol, [s_x(l).xi], [s_y(l).xi]));
        case {'clade', 'unclade', 'parent', 'sibling', 'child'}
            ts(i) = isequaln({s_x.(fi)}, {s_y.(fi)});
        case 'time'
            ts(i) = all(arrayfun(@ismembertol, [s_x([l, n]).time], ...
                                 [s_y([l, n]).time]));
        case 'catloc'
            ts(i) = all(cellfun(@(x, y) all(arrayfun(@ismembertol, x, y)), ...
                                {s_x(state_x.cat > 0).catloc}, ...
                                {s_y(state_y.cat > 0).catloc}));
        otherwise
            tv(i) = 0;
        end

    end
    t = all(ts(tv == 1));
end
