function [state_x, succ_x, state_y, succ_y] = MarkovCoupledMaximal(mcmc, ...
        model, state_x, state_y, ignoreearlywarn, MV, u_mh)
    % Maximal coupling of scalar parameter proposal and acceptance steps

    global NARROW WIDE DEPNU VARYMU DONTMOVECATS BORROWING VARYBETA

    OK_x = 1;
    OK_y = 1;
    nstate_x = state_x;
    nstate_y = state_y;

    TOPOLOGY = 0;

    switch MV
        case 1
            update='RW node time between parent time and oldest child time';
            [i, newage_x, newage_y, logq_x, logq_y] ...
                = SchooseCoupledMaximal(state_x, state_y);
            % Supdate always returns TOPOLOGY = 0 so we ignore it
            [nstate_x, U_x, ~] = Supdate(state_x, i, newage_x);
            [nstate_y, U_y, ~] = Supdate(state_y, i, newage_y);

            if BORROWING
                logq_x = logq_x + catastropheScalingFactor(state_x, nstate_x);
                logq_y = logq_y + catastropheScalingFactor(state_y, nstate_y);
            end
        case {4, 5}
            if MV == 4
                update = 'Reconnect an edge into a nearby edge';
                mt = NARROW;
            else
                update = 'Reconnect an edge into an edge chosen UAR over the tree';
                mt = WIDE;
            end
            [i, j_x, j_y, k_x, k_y, newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupledMaximal(state_x, state_y, mt, ...
                                        mcmc.update.theta, model.prior);
            % Bupdate always returns TOPOLOGY = 1 so we ignore it below
            TOPOLOGY = 1;
            OK_x = ~isempty(newage_x);
            if OK_x
                [nstate_x, U_x, ~] = Bupdate(state_x, i, j_x, k_x, newage_x);
                if BORROWING
                    logq_x = logq_x ...
                             + catastropheScalingFactor(state_x, nstate_x);
                end
            end
            OK_y = ~isempty(newage_y);
            if OK_y
                [nstate_y, U_y, ~] = Bupdate(state_y, i, j_y, k_y, newage_y);
                if BORROWING
                    logq_y = logq_y ...
                             + catastropheScalingFactor(state_y, nstate_y);
                end
            end
        case 8
            update = 'Vary mu';
            if ~VARYMU
                disp('vary mu?');
                keyboard;
            end

            [nstate_x.mu, nstate_y.mu] ...
                = maximalCouplingUniformScaling(...
                    state_x.mu, ...
                    state_y.mu, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

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

            [nstate_x.rho, nstate_y.rho] ...
                = maximalCouplingUniformScaling(...
                    state_x.rho, ...
                    state_y.rho, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

            var_x = nstate_x.rho / state_x.rho;
            var_y = nstate_y.rho / state_y.rho;

            logq_x = -log(var_x);
            logq_y = -log(var_y);

            U_x = [];
            U_y = [];
        case 16
            update = 'Vary kappa';

            [nstate_x.kappa, nstate_y.kappa] ...
                = maximalCouplingUniformScaling(...
                    state_x.kappa, ...
                    state_y.kappa, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

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

            [nstate_x.lambda, nstate_y.lambda] ...
                = maximalCouplingUniformScaling(...
                    state_x.lambda, ...
                    state_y.lambda, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

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
            leaf = randsample(state_x.leaves, 1);
            % state_x.leaves(ceil(rand * length(state.leaves)));

            % xi'_tilde = 1 - xi' = U(del, del + deldel) * (1 - xi)
            [nstate_x_tree_leaf_xi_tilde, nstate_y_tree_leaf_xi_tilde] ...
                = maximalCouplingUniformScaling(...
                    1 - state_x.tree(leaf).xi, ...
                    1 - state_y.tree(leaf).xi, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

            nstate_x.tree(leaf).xi = 1 - nstate_x_tree_leaf_xi_tilde;
            nstate_y.tree(leaf).xi = 1 - nstate_y_tree_leaf_xi_tilde;

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

            [nstate_x.beta, nstate_y.beta] ...
                = maximalCouplingUniformScaling(...
                    state_x.beta, ...
                    state_y.beta, ...
                    mcmc.update.del, ...
                    mcmc.update.del + mcmc.update.deldel);

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
