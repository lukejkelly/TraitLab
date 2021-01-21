function [i, j_x, j_y, k_x, k_y, newage_x, newage_y, logq_x, logq_y] ...
        = BchooseCoupledMaximal(state_x, state_y, mt, THETA, prior)

    % Housekeeping means that Adam, root, clade roots and leaves have same
    % indices, internal nodes within a clade have a common set of indices

    global ROOT NARROW WIDE

    N = 2 * state_x.NS - 1;
    s_x = state_x.tree;
    s_y = state_y.tree;

    % Part I: picking same edge to move in both trees
    i = ceil(N * rand);
    while s_x(i).type == ROOT
       i = ceil(N * rand);
    end

    % Part II: picking destination
    switch mt
    case NARROW
        % Narrow, so move <iP, i> into <pa(iP), sib(iP)> = <k, j>
        [j_x, k_x, FAIL_x] = BchooseCoupledMaximal.getNarrowDestination(i, s_x);
        [j_y, k_y, FAIL_y] = BchooseCoupledMaximal.getNarrowDestination(i, s_y);
    case WIDE
        % Choose destination <k, j> uniformly across tree / clade iP
        if prior.isclade
            r_x = BchooseCoupledMaximal.getWideCandidatesClade(i, s_x);
            r_y = BchooseCoupledMaximal.getWideCandidatesClade(i, s_y);
            if ~isequaln(r_x, r_y)
                disp(r_x);
                disp(r_y);
                save(sprintf('output/%d.mat', yyyymmdd(datetime)));
                error('Destination edge sets should be identical');
            end
            r = r_x;
            N = length(r);
        else
            r = 1:N; % TODO: Why does this allow root but above does not?
        end
        r = r(randperm(N));
        [j_x, k_x, FAIL_x] ...
            = BchooseCoupledMaximal.getWideDestination(i, r, N, s_x);
        [j_y, k_y, FAIL_y] ...
            = BchooseCoupledMaximal.getWideDestination(i, r, N, s_y);

        % If j_x is the root in x then j_y = j_x is the root in y
        % TODO: remove/correct this check as j is only the root if it appears
        % before any other valid option in r
        % if ~(FAIL_x || FAIL_y) ...
        %         && any([s_x(j_x).type, s_y(j_y).type] == ROOT) ...
        %         && (j_x ~= j_y || s_x(j_x).type ~= s_y(j_y).type)
        %     disp([j_x, j_y]);
        %     error('Root indices should match');
        % end
    otherwise
        error('Move type must be NARROW or WIDE')
    end

    % Part III: sample new node time
    if ~(FAIL_x || FAIL_y)
        % Attempt to couple node times
        [newage_x, newage_y, logq_x, logq_y] ...
            = BchooseCoupledMaximal.sampleCoupling(i, j_x, j_y, k_x, k_y, ...
                                                   s_x, s_y, THETA);
    else
        % Attempt to sample from marginal node times
        if FAIL_x
            newage_x = [];
            logq_x = -Inf;
        else
            [newage_x, logq_x] ...
                = BchooseCoupledMaximal.sampleMarginal(i, j_x, k_x, s_x, THETA);
        end
        if FAIL_y
            newage_y = [];
            logq_y = -Inf;
        else
            [newage_y, logq_y] ...
                = BchooseCoupledMaximal.sampleMarginal(i, j_y, k_y, s_y, THETA);
        end
    end
end
