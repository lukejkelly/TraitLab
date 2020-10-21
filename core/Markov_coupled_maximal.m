function [state_x, succ_x, state_y, succ_y] = Markov_coupled_maximal(mcmc, ...
    model, state_x, state_y, ignoreearlywarn, MV, u_mh)
    % Maximal coupling of scalar parameter proposal and acceptance steps
    % Currently MV = 8 (mu) or 21 (beta)

    global DEPNU VARYMU DONTMOVECATS BORROWING VARYBETA

    OK_x = 1;
    OK_y = 1;
    nstate_x = state_x;
    nstate_y = state_y;

    TOPOLOGY = 0;

    % Prototypes for coupling
    rp = @(v) v * (mcmc.update.del + rand * mcmc.update.deldel);
    dp = @(u, v) (v * mcmc.update.del <= u ...
                  && u <= v * (mcmc.update.del + mcmc.update.deldel)) ...
                 / (v * mcmc.update.deldel);

    switch MV
        case 8
            update='Vary mu';
            if ~VARYMU
                disp('vary mu?');
                keyboard;
            end

            rx = @() rp(state_x.mu);
            dx = @(x) dp(x, state_x.mu);
            ry = @() rp(state_y.mu);
            dy = @(y) dp(y, state_y.mu);

            [nstate_x.mu, nstate_y.mu] = coupling_maximal(rx, dx, ry, dy);

            var_x = nstate_x.mu / state_x.mu;
            var_y = nstate_y.mu / state_y.mu;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = state_x.nodes;
            U_y = state_y.nodes;

            if DEPNU
                nstate_x.nu = state_x.nu / var_x;
                nstate_y.nu = state_y.nu / var_y;
            end
            if DONTMOVECATS && nstate_x.mu < 1e-5
                OK_x = 0;
            end
            if DONTMOVECATS && nstate_y.mu < 1e-5
                OK_y = 0;
            end
        case 21
            update = 'Vary beta';
            if ~VARYBETA && BORROWING
                disp('vary beta?')
                keyboard;
            end

            rx = @() rp(state_x.beta);
            dx = @(x) dp(x, state_x.beta);
            ry = @() rp(state_y.beta);
            dy = @(y) dp(y, state_y.beta);

            [nstate_x.beta, nstate_y.beta] = coupling_maximal(rx, dx, ry, dy);

            var_x = nstate_x.beta / state_x.beta;
            var_y = nstate_y.beta / state_y.beta;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = state_x.nodes;
            U_y = state_y.nodes;
    end

    [state_x, succ_x] = updateState(update, model, u_mh, TOPOLOGY, logq_x, ...
        state_x, nstate_x, OK_x, U_x, ignoreearlywarn);
    [state_y, succ_y] = updateState(update, model, u_mh, TOPOLOGY, logq_y, ...
        state_y, nstate_y, OK_y, U_y, ignoreearlywarn);
end

function [state, succ] = updateState(update, model, u_mh, TOPOLOGY, ...
    logq, state, nstate, OK, U, ignoreearlywarn)

    global TESTUP BORROWING

    if OK && model.prior.isclade
        if TOPOLOGY
            nstate = UpdateClades(nstate, U, size(model.prior.clade, 2));
        end
        OK = CladePrior(model.prior, nstate);
    end

    succ = 0; % Switched to 1 if move accepted

    if OK
        if ~all(sort(find([state.tree.type] == 0)) == sort(find([nstate.tree.type] == 0)))
            disp('Leaf variables mismatch')
            keyboard;
        end

        if TOPOLOGY && DONTMOVECATS
            nstate.cat(:) = 1;
            nstate.cat(nstate.root) = 0;
            nstate.cat(nstate.tree(nstate.root).parent) = 0;
            nstate.ncat = sum(nstate.cat);
        end
        [nstate, ~] = MarkRcurs(nstate, U, TOPOLOGY, ignoreearlywarn);

        % Likelihood calculations.
        if BORROWING
            [nstate.loglkd, nstate.fullloglkd] = logLkd2(nstate);
        else
            nstate.loglkd = LogLkd(nstate);
            nstate.fullloglkd = LogLkd(nstate, nstate.lambda);
        end

        % Log-prior for tree and catastrophes
        nstate.logprior = LogPrior(model.prior, nstate);

        % Log-prior for parameters: mu, beta, etc.
        logpp = LogPriorParm(state, nstate);

        % Log-Hastings ratio in acceptance step
        logr = nstate.logprior - state.logprior + nstate.loglkd - state.loglkd + logq + logpp;

        % Acceptance step
        if ( (logr > 0) || (log(u_mh) < logr) || isinf(state.loglkd) )
            state = nstate;
            succ = 1;
        end

        if TESTUP && check(state, [])
            disp(['Error in update:', update]);
            check(nstate, [])
            keyboard;
        end
    end
end

% Sampling from a maximal coupling of two distributions
function [x, y] = coupling_maximal(rp, dp, rq, dq)
    x = rp();
    if rand * dp(x) <= dq(x)
        y = x;
    else
        y = nan;
        while isnan(y)
            ys = rq();
            if rand * dq(ys) > dp(ys)
                y = ys;
            end
        end
    end
end
