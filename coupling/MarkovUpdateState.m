function [state, succ] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
        logq, state, nstate, OK, U, ignoreearlywarn)

    global TESTUP BORROWING DONTMOVECATS

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
