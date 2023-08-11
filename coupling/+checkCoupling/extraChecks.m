function t = extraChecks(state_x, state_y)
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
        case {'leaves', 'nodes', 'root'}
            ts(i) = all(state_x.(fi) == state_y.(fi));
            if ts(i) == 0
                error('Housekeeping failure on %s', fi);
            end
        case 'claderoot'
            ts(i) = all(isequaln(state_x.(fi), state_y.(fi)));
            if ts(i) == 0
                error('Housekeeping failure on %s', fi);
            end
        case {'mu', 'beta', 'kappa'} % 'lambda', 'rho'
            ts(i) = ismembertol(state_x.(fi), state_y.(fi));
        case {'loglkd', 'logprior'}
            ts(i) = ismembertol(state_x.(fi), state_y.(fi));
        case 'tree'
            ts(i) = compareTrees(state_x, state_y);
        case 'length'
            ts(i) = ismembertol(state_x.(fi), state_y.(fi));
        case {'cat', 'ncat'}
            ts(i) = all(state_x.(fi) == state_y.(fi));
        otherwise
            % Do not care for p, fullloglkd or nu
            ts(i) = 1;
        end
        if ~ts(i)
            fprintf('state.%s not matching\n', fi);
            printNotMatching('x', fi, state_x);
            printNotMatching('y', fi, state_y);
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
        case {'u', 'v', 'LamInt', 'w', 'PD', 'ActI', 'CovI', 'difCovI', ...
              'Tu', 'TActI', 'TCovI', 'TdifCovI', 'mark', ...
              'PDcat', 'CatLamInt', 'n', 'd'}
            % Used to calculate likelihood in SD so can be ignored
            % t(i) = all(cellfun(@isequaln, {s_x.(fi)}, {s_y.(fi)}));
            ts(i) = 1;
        case {'cat', 'catloc'}
            % Handled differently by SD and SDLT so ignoring for now
            ts(i) = 1;
        case {'time'}
            ts(i) = all(cellfun(@ismembertol, {s_x.(fi)}, {s_y.(fi)}));
        otherwise
            ts(i) = isequaln({s_x.(fi)}, {s_y.(fi)});
        end
        if ~ts(i)
            fprintf('tree.%s not matching\n', fi);
            printNotMatching('x', fi, s_x);
            printNotMatching('y', fi, s_y);
        end
    end
    t = all(ts);
end

function printNotMatching(d, f, s)
    t = formattedDisplayText({s.(f)}, 'NumericFormat', 'longg');
    fprintf('_%s.%s\n%s\n', d, f, t);
end
