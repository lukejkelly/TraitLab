function [state_x, succ_x, state_y, succ_y] = MarkovCoupledMaximal(mcmc, ...
    model, state_x, state_y, ignoreearlywarn, MV, u_mh)
    % Maximal coupling of scalar parameter proposal and acceptance steps

    global DEPNU VARYMU DONTMOVECATS BORROWING VARYBETA

    OK_x = 1;
    OK_y = 1;
    nstate_x = state_x;
    nstate_y = state_y;

    TOPOLOGY = 0;

    % Prototypes for coupling proposals
    rp = @(v) v * (mcmc.update.del + rand * mcmc.update.deldel);
    dp = @(u, v) (v * mcmc.update.del <= u ...
                  && u <= v * (mcmc.update.del + mcmc.update.deldel)) ...
                 / (v * mcmc.update.deldel);
    switch MV
        case 8
            update = 'Vary mu';
            if ~VARYMU
                disp('vary mu?');
                keyboard;
            end

            rx = @() rp(state_x.mu);
            dx = @(x) dp(x, state_x.mu);
            ry = @() rp(state_y.mu);
            dy = @(y) dp(y, state_y.mu);

            [nstate_x.mu, nstate_y.mu] = maximalCoupling(rx, dx, ry, dy);

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

            OK_x = ~(DONTMOVECATS && nstate_x.mu < 1e-5);
            OK_y = ~(DONTMOVECATS && nstate_y.mu < 1e-5);
        case 15
            update = 'Vary rho';

            rx = @() rp(state_x.rho);
            dx = @(x) dp(x, state_x.rho);
            ry = @() rp(state_y.rho);
            dy = @(y) dp(y, state_y.rho);

            [nstate_x.rho, nstate_y.rho] = maximalCoupling(rx, dx, ry, dy);

            var_x = nstate_x.rho / state_x.rho;
            var_y = nstate_y.rho / state_y.rho;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = [];
            U_y = [];
        case 16
            update = 'Vary kappa';

            rx = @() rp(state_x.kappa);
            dx = @(x) dp(x, state_x.kappa);
            ry = @() rp(state_y.kappa);
            dy = @(y) dp(y, state_y.kappa);

            [nstate_x.kappa, nstate_y.kappa] = maximalCoupling(rx, dx, ry, dy);

            var_x = nstate_x.kappa / state_x.kappa;
            var_y = nstate_y.kappa / state_y.kappa;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = state_x.nodes;
            U_y = state_y.nodes;

            OK_x = ((0.25 <= nstate_x.kappa) && (nstate_x.kappa <= 1));
            OK_y = ((0.25 <= nstate_y.kappa) && (nstate_y.kappa <= 1));

            if OK_x && DEPNU
                nstate_x.nu = state_x.nu / var_x;
            end
            if OK_y && DEPNU
                nstate_y.nu = state_y.nu / var_y;
            end
        case 17
            update = 'Vary lambda';

            rx = @() rp(state_x.lambda);
            dx = @(x) dp(x, state_x.lambda);
            ry = @() rp(state_y.lambda);
            dy = @(y) dp(y, state_y.lambda);

            [nstate_x.lambda, nstate_y.lambda] = maximalCoupling(rx, dx, ry, dy);

            var_x = nstate_x.lambda / state_x.lambda;
            var_y = nstate_y.lambda / state_y.lambda;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = state_x.nodes;
            U_y = state_y.nodes;

            if DEPNU
                nstate_x.nu=state_x.nu * var_x;
                nstate_y.nu=state_y.nu * var_y;
            end
        case 19
            update = 'Vary XI for one leaf';
            % Same leaf indices in both x and y after housekeeping
            leaf = randsample(state_x.leaves, 1); % state_x.leaves(ceil(rand * length(state.leaves)));

            rx = @() 1 - rp(1 - state_x.tree(leaf).xi);
            dx = @(x) dp(1 - x, 1 - state_x.tree(leaf).xi);
            ry = @() 1 - rp(1 - state_y.tree(leaf).xi);
            dy = @(y) dp(1 - y, 1 - state_y.tree(leaf).xi);

            [nstate_x.tree(leaf).xi, nstate_y.tree(leaf).xi] ...
                = maximalCoupling(rx, dx, ry, dy);

            var_x = (1 - nstate_x.tree(leaf).xi) / (1 - state_x.tree(leaf).xi);
            var_y = (1 - nstate_y.tree(leaf).xi) / (1 - state_y.tree(leaf).xi);

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = above(leaf, state_x.tree, state_x.root);
            U_y = above(leaf, state_y.tree, state_y.root);

            OK_x = nstate_x.tree(leaf).xi >= 0;
            OK_y = nstate_y.tree(leaf).xi >= 0;
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

            [nstate_x.beta, nstate_y.beta] = maximalCoupling(rx, dx, ry, dy);

            var_x = nstate_x.beta / state_x.beta;
            var_y = nstate_y.beta / state_y.beta;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = state_x.nodes;
            U_y = state_y.nodes;
        otherwise
            error('Not implemented');
    end

    [state_x, succ_x] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
        logq_x, state_x, nstate_x, OK_x, U_x, ignoreearlywarn);
    [state_y, succ_y] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
        logq_y, state_y, nstate_y, OK_y, U_y, ignoreearlywarn);
end
