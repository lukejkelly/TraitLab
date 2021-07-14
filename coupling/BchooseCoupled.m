function [i, j_x, j_y, k_x, k_y, newage_x, newage_y, logq_x, logq_y] ...
        = BchooseCoupled(state_x, state_y, mt, THETA, prior)

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
        [j_x, k_x, FAIL_x] = BchooseCoupled.getNarrowDestination(i, s_x);
        [j_y, k_y, FAIL_y] = BchooseCoupled.getNarrowDestination(i, s_y);
    case WIDE
        % Choose destination <k, j> uniformly across tree / clade iP
        if prior.isclade
            r_x = BchooseCoupled.getWideCandidatesClade(i, s_x);
            r_y = BchooseCoupled.getWideCandidatesClade(i, s_y);
            if ~isequaln(r_x, r_y)
                disp(r_x);
                disp(r_y);
                save(sprintf('output/%d.mat', yyyymmdd(datetime)));
                error('Destination edge sets should be identical');
            end
            r = r_x;
            N = length(r);
        else
            r = 1:N;
        end

        if  N > 4
            [j_x, j_y, k_x, k_y, FAIL_x, FAIL_y] ...
                = BchooseCoupled.getWideDestination(i, r, s_x, s_y);
        else
            [j_x, k_x, FAIL_x] = deal(-1, -1, 1);
            [j_y, k_y, FAIL_y] = deal(-1, -1, 1);
        end
    otherwise
        error('Move type must be NARROW or WIDE')
    end

    % Part III: sample new node time
    if ~(FAIL_x || FAIL_y)
        % Attempt to couple node times
        [newage_x, newage_y, logq_x, logq_y] ...
            = BchooseCoupled.sampleCoupling(i, j_x, j_y, k_x, k_y, s_x, s_y, ...
                                            THETA);
    else
        % Attempt to sample from marginal node times
        if FAIL_x
            newage_x = [];
            logq_x = -Inf;
        else
            [newage_x, logq_x] ...
                = BchooseCoupled.sampleMarginal(i, j_x, k_x, s_x, THETA);
        end
        if FAIL_y
            newage_y = [];
            logq_y = -Inf;
        else
            [newage_y, logq_y] ...
                = BchooseCoupled.sampleMarginal(i, j_y, k_y, s_y, THETA);
        end
    end
end
