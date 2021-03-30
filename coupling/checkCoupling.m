function coupled = checkCoupling(state_x, state_y)
    % TODO: Check catastrophes
    coupled = state_x.logprior == state_y.logprior ...
        && state_x.beta == state_y.beta ...
        && state_x.loglkd == state_y.loglkd ...
        && state_x.mu == state_y.mu ...
        && all([state_x.tree.time] == [state_y.tree.time]) ...
        && equaltrees(state_x.tree, state_y.tree);

    % TODO: change after further checks
    if coupled && ~compareStates(state_x, state_y)
        warning('Some state fields do not match');
        save(sprintf('coupling/compare-state %s.mat', datetime));
    end
end

function t = compareStates(state_x, state_y)
    % Extra state checks
    fn = fieldnames(state_x);
    ts = zeros(size(fn));
    for i = 1:length(ts)
        fi = fn{i};
        switch fi
        case 'tree'
            ts(i) = compareTrees(state_x, state_y);
        otherwise
            ts(i) = isequaln({state_x.(fi)}, {state_y.(fi)});
        end
    end
    t = all(ts);
end

function t = compareTrees(state_x, state_y)
    % Extra tree checks
    [s_x, l_x] = deal(state_x.tree, state_x.leaves);
    [s_y, l_y] = deal(state_y.tree, state_y.leaves);
    fn = fieldnames(s_x);
    ts = zeros(size(fn));
    for i = 1:length(ts)
        fi = fn{i};
        switch fi
        case 'Name'
            ts(i) = all(strcmp({s_x(l_x).(fi)}, {s_y(l_y).(fi)}));
        case 'xi'
            ts(i) = isequaln({s_x(l_x).(fi)}, {s_y(l_x).(fi)});
        case {'clade', 'unclade'}
            ts(i) = all(cellfun(@(x, y) (isempty(x) && isempty(y)) ...
                                        || isequaln(x, y), ...
                                {s_x.(fi)}, {s_y.(fi)}));
        case {'w', 'PD', 'TCovI', 'PDcat'}
            % Used to calculate likelihood in SD so can be ignored
            % t(i) = all(cellfun(@isequaln, {s_x.(fi)}, {s_y.(fi)}));
            ts(i) = 1;
        otherwise
            ts(i) = isequaln({s_x.(fi)}, {s_y.(fi)});
        end
    end
    if any(~ts)
        disp(fn(~ts));
    end
    t = all(ts);
end
