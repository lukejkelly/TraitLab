function [nstate_x, nstate_y, logq_x, logq_y, U_x, U_y] ...
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
        [j_x, k_x, FAIL_x] = getNarrowDestination(i, s_x);
        [j_y, k_y, FAIL_y] = getNarrowDestination(i, s_y);
    case WIDE
        % Choose destination <k, j> uniformly across tree / clade iP
        if prior.isclade
            r_x = getWideCandidatesClade(i, s_x);
            r_y = getWideCandidatesClade(i, s_y);
            if ~isequaln(r_x, r_y)
                disp(r_x);
                disp(r_y);
                error('Destination edge sets should be identical');
            end
            r = r_x;
            N = length(r);
        else
            r = 1:N;
        end
        r = r(randperm(N));
        [j_x, k_x, FAIL_x] = getWideDestination(i, r, N, s_x);
        [j_y, k_y, FAIL_y] = getWideDestination(i, r, N, s_y);

        % If j_x is the root in x then j_y = j_x is the root in y
        if any([s_x(j_x).type, s_x(j_x).type] == ROOT) ...
            && (j_x ~= j_y || s_x(j_x).type ~= s_x(j_x).type)
            disp([j_x, j_y]);
            error('Root indices should match');
        end
    end

    % Part III: sample new node time
    if ~(FAIL_x || FAIL_y)
        % Attempt to couple node times
        [newage_x, newage_y, logq_x, logq_y] ...
            = sampleCoupling(i, j_x, j_y, k_x, k_y, s_x, s_y, THETA);
    else
        % Attempt to sample from marginal node times
        if ~FAIL_x
            [newage_x, logq_x] = sampleMarginal(i, j_x, k_x, s_x, THETA);
        end
        if ~FAIL_y
            [newage_y, logq_y] = sampleMarginal(i, j_y, k_y, s_y, THETA);
        end
    end

    % Part IV: create new state
    if FAIL_x
        nstate_x = [];
        logq_x = -Inf;
        U_x = [];
    else
        [nstate_x, U_x] = Bupdate(state_x, i, j_x, k_x, newage_x);
    end
    if FAIL_y
        nstate_y = [];
        logq_y = [];
        U_y = [];
    else
        [nstate_y, U_y] = Bupdate(state_y, i, j_y, k_y, newage_y);
    end
end

function [j, k, FAIL] = getNarrowDestination(i, s)
    % <k, j> = <pa(pa(i)), sib(pa(i))>
    global ROOT OTHER
    iP = s(i).parent;
    if s(iP).type == ROOT
        FAIL = 1;
        k = -1;
        j = -1;
    else
        FAIL = 0;
        k = s(iP).parent;
        j = s(k).child(OTHER(s(iP).sibling));
    end
end

function r = getWideCandidatesClade(i, s)
    % Possible destination edges for a wide SPR move
    global ROOT
    r = [];
    iP = s(i).parent;
    for a = 1:length(s)
        % I think what this is checking is whether cl(pa(a)) is equal to
        % cl(pa(i)) which is a subset (not strict) of cl(a):
        %   |cl(a)| >= |cl(pa(i))| AND cl(pa(i)) \ cl(a) = EMPTYSET
        %       -> cl(pa(i)) SUBSET cl(a)
        %   |cl(pa(a))| >= |cl(pa(i))| AND cl(pa(a)) \ cl(pa(i)) = EMPTYSET
        %       -> cl(pa(i)) = cl(a)
        if s(a).type < ROOT ...
                && length(s(a).clade) >= length(s(iP).clade) ...
                && length(s(s(a).parent).clade) >= length(s(iP).clade) ...
                && isempty(setdiff(s(iP).clade, s(a).clade)) ...
                && isempty(setdiff(s(s(a).parent).clade, s(iP).clade))
            r = [r, a];
        end
    end
    % TODO: can the above be replaced by an array function
    r2 = find(arrayfun(@(a) s(a).type < ROOT ...
                       && length(s(a).clade) >= length(s(iP).clade) ...
                       && length(s(s(a).parent).clade) >= length(s(iP).clade) ...
                       && isempty(setdiff(s(iP).clade, s(a).clade)) ...
                       && isempty(setdiff(s(s(a).parent).clade, s(iP).clade)), ...
                       1:length(s)));
    if ~isequaln(r, r2)
        warning('Fix array function');
    end
    % TODO: can we simplify into set equality and inclusion?
    r3 = find(arrayfun(@(a) s(a).type < ROOT ...
                       && all(ismember(s(iP).clade, s(a).clade)) ...
                       && isequaln(s(s(a).parent).clade, s(iP).clade), ...
                       1:length(s)));
    if ~isequaln(r, r3)
        warning('Fix array function');
    end
end

function [j, k, FAIL] = getWideDestination(i, r, N, s)
    if N > 4
        iT = s(i).time;
        rInd = 1;
        j = r(rInd);
        k = s(j).parent;
        while (s(k).time <= iT || i == j || i == k)
            rInd = rInd + 1;
            j = r(rInd);
            k = s(j).parent;
        end
        FAIL = any([j, k] == s(i).parent);
    else
        j = -1;
        k = -1;
        FAIL = 1;
    end
end

function [newage, logq] = sampleMarginal(i, j, k, s, THETA)
    global OTHER ROOT

    iT = s(i).time;
    iP = s(i).parent;

    if s(j).type == ROOT

       jT = s(j).time;
       delta = -(1 / THETA) * log(rand);
       newage = jT + delta;

       PiP = s(iP).parent;
       CiP = s(iP).child(OTHER(s(i).sibling));

       PiPT = s(PiP).time;
       CiPT = s(CiP).time;
       old_minage = max(iT,CiPT);
       old_range = PiPT - old_minage;

       q = exp(delta * THETA) / (THETA * old_range);
       logq = log(q);

    elseif s(iP).type == ROOT

       jT = s(j).time;
       kT = s(k).time;
       new_minage = max(iT, jT);
       new_range = kT - new_minage;
       newage = new_minage + rand * new_range;

       CiP = s(iP).child(OTHER(s(i).sibling));
       CiPT = s(CiP).time;
       q = exp((CiPT - s(iP).time) * THETA) * new_range * THETA;
       logq = log(q);

    else

       jT = s(j).time;
       kT = s(k).time;
       new_minage = max(iT, jT);
       new_range = kT - new_minage;
       newage = new_minage + rand * new_range;

       PiP = s(iP).parent;
       CiP = s(iP).child(OTHER(s(i).sibling));

       PiPT = s(PiP).time;
       CiPT = s(CiP).time;
       old_minage = max(iT, CiPT);
       old_range = PiPT - old_minage;

       q = new_range / old_range;
       logq = log(q);
    end
end

function [newage_x, newage_y, logq_x, logq_y] = sampleCoupling(i, j_x, j_y, ...
        k_x, k_y, s_x, s_y, THETA)

    global ROOT

    if s_x(j_x).type == ROOT
        % If j_x is the root in x then j_y = j_x and is the root in y
        jT_x = s_x(j_x).time;
        jT_y = s_y(j_y).time;
        [newage_x, newage_y] ...
            = maximalCouplingShiftedExponential(jT_x, jT_y, THETA);
        logq_x = getCoupling1Parameter(i, jT_x, s_x, newage_x, THETA);
        logq_y = getCoupling1Parameter(i, jT_y, s_y, newage_y, THETA);
    else
        % logq depends on whether iP is the root or not
        [new_minage_x, kT_x, logq_x] ...
            = getCoupling2Parameters(i, j_x, k_x, s_x, THETA);
        [new_minage_y, kT_y, logq_y] ...
            = getCoupling2Parameters(i, j_y, k_y, s_y, THETA);
        [newage_x, newage_y] ...
            = maximalCouplingUniform(new_minage_x, kT_x, new_minage_y, kT_y);
    end
end

function [logq] = getCoupling1Parameter(i, jT, s, newage, THETA)
    global OTHER

    delta = newage - jT;

    iT = s(i).time;
    iP = s(i).parent;

    PiP = s(iP).parent;
    CiP = s(iP).child(OTHER(s(i).sibling));

    PiPT = s(PiP).time;
    CiPT = s(CiP).time;

    old_minage = max(iT, CiPT);
    old_range = PiPT - old_minage;

    logq = delta * THETA - log(THETA * old_range);
end

function [new_minage, kT, logq] = getCoupling2Parameters(i, j, k, s, THETA)
    global OTHER ROOT

    iT = s(i).time;
    jT = s(j).time;
    kT = s(k).time;

    iP = s(i).parent;

    new_minage = max(iT, jT);
    new_range = kT - new_minage;

    PiP = s(iP).parent;
    CiP = s(iP).child(OTHER(s(i).sibling));

    PiPT = s(PiP).time;
    CiPT = s(CiP).time;
    old_minage = max(iT, CiPT);
    old_range = PiPT - old_minage;

    % Sampling is the same but logq depends on whether iP was root
    if s(PiP).type == ROOT
        logq = (CiPT - s(iP).time) * THETA + log(new_range) + log(THETA);
    else
        logq = log(new_range) - log(old_range);
    end
end
