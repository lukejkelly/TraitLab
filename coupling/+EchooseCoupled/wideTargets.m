function v = wideTargets(state, prior)
    % Finding valid combinations of i and j, based on core/Echoose
    global ROOT OTHER

    s = state.tree;
    Root = state.root;
    N = 2 * state.NS - 1;

    u = nan(N);
    for i = 1:N
        if (i == Root || (s(s(i).parent).type == ROOT && s(s(s(i).parent).child(OTHER(s(i).sibling))).time < s(i).time))
            % i not valid if root or root's eldest child
            u(i, :) = 0;
        else
            iP = s(i).parent;
            for j = 1:N
                jP = s(j).parent;
                if i ~= j && iP ~= jP && i ~= jP && j ~= iP && s(j).time < s(iP).time && s(i).time < s(jP).time ...
                        && (~prior.isclade || (prior.isclade && isequal(s(iP).clade, s(jP).clade)))
                    u(i, j) = 1;
                else
                    u(i, j) = 0;
                end
            end
        end
    end
    v = find(u);
end
