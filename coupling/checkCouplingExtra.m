function t = checkCouplingExtra(state_x, state_y)
    % Extra state checks
    fn = fieldnames(state_x);
    ts = zeros(size(fn));
    for i = 1:length(ts)
        fi = fn{i};
        switch fi
        case {'NS', 'L'}
            if state_x.(fi) == state_y.(fi)
                ts(i) = 1;
            else
                error('Leaf and data counts do not match');
            end
        case {'leaves', 'root'}
            if all(sort(state_x.(fi)) == sort(state_y.(fi)))
                ts(i) = 1;
            else
                error('Housekeeping failure');
            end
        case 'tree'
            ts(i) = compareTrees(state_x, state_y);
        otherwise
            ts(i) = isequaln({state_x.(fi)}, {state_y.(fi)});
        end
    end
    if any(~ts)
        fprintf('state fields not matching:\t%s\n', strjoin(fn(~ts)));
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
        case {'u', 'v', 'LamInt', 'w', 'PD', 'ActI', 'CovI', 'difCovI', ...
              'Tu', 'TActI', 'TCovI', 'TdifCovI', 'mark', ...
              'PDcat', 'CatLamInt', 'n', 'd'}
            % Used to calculate likelihood in SD so can be ignored
            % t(i) = all(cellfun(@isequaln, {s_x.(fi)}, {s_y.(fi)}));
            ts(i) = 1;
        case {'cat', 'catloc'}
            % Handled differently by SD and SDLT so ignoring for now
            ts(i) = 1;
        otherwise
            ts(i) = isequaln({s_x.(fi)}, {s_y.(fi)});
        end
    end
    if any(~ts)
        fprintf('tree fields not matching:\t%s\n', strjoin(fn(~ts)));
    end
    t = all(ts);
end
