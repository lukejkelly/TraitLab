function [state_x, succ_x, state_y, succ_y] = MarkovCoupledMaximal(mcmc, ...
        model, state_x, state_y, ignoreearlywarn, MV, u_mh)
    % Maximal coupling of scalar parameter proposal and acceptance steps

    global NARROW WIDE VARYMU BORROWING VARYBETA

    OK_x = 1;
    OK_y = 1;
    TOPOLOGY = 0;

    switch MV
    case 1
        update='RW node time between parent time and oldest child time';
        [i, newage_x, newage_y, logq_x, logq_y] ...
            = SchooseCoupled(state_x, state_y);
        % Supdate always returns TOPOLOGY = 0 so we ignore it
        [nstate_x, U_x, ~] = Supdate(state_x, i, newage_x);
        [nstate_y, U_y, ~] = Supdate(state_y, i, newage_y);
    case {4, 5}
        if MV == 4
            update = 'Reconnect an edge into a nearby edge';
            mt = NARROW;
        else
            update = 'Reconnect an edge into an edge chosen UAR over the tree';
            mt = WIDE;
        end
        [i, j_x, j_y, k_x, k_y, newage_x, newage_y, logq_x, logq_y] ...
            = BchooseCoupled(state_x, state_y, mt, mcmc.update.theta, ...
                             model.prior);
        % Bupdate always returns TOPOLOGY = 1 so we ignore it below
        TOPOLOGY = 1;
        OK_x = ~isempty(newage_x);
        if OK_x
            [nstate_x, U_x, ~] = Bupdate(state_x, i, j_x, k_x, newage_x);
        end
        OK_y = ~isempty(newage_y);
        if OK_y
            [nstate_y, U_y, ~] = Bupdate(state_y, i, j_y, k_y, newage_y);
        end
    case 6
        update = 'Rescale whole tree trying to match root times';
        [var_x, var_y] = RscaleCoupled(state_x, state_y, mcmc.update.del, ...
                                       mcmc.update.del + mcmc.update.deldel);
        [nstate_x, U_x, ~, OK_x, logq_x] = Rscale(state_x, var_x);
        [nstate_y, U_y, ~, OK_y, logq_y] = Rscale(state_y, var_y);
    case 7
        update='Rescale randomly chosen subtree';
        [nstate_x, nstate_y, U_x, U_y, logq_x, logq_y, OK_x, OK_y] ...
            = RscaleSubTreeCoupled(state_x, state_y, mcmc.update.del, ...
                                   mcmc.update.deldel);
    case 8
        update = 'Vary mu';
        if ~VARYMU
            disp('vary mu?');
            keyboard;
        end
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = scaleCoupledMu(state_x, state_y, mcmc);
    case 12
        update = 'Rescale top tree';
        if ~model.prior.isclade
            error('Move 12 should not be selected if no clades');
        end
        [nstate_x, nstate_y, U_x, U_y, logq_x, logq_y, OK_x, OK_y] ...
            = RscaleTopTreeCoupled(state_x, state_y, model.prior, ...
                                   mcmc.update.del, mcmc.update.deldel);
    case 13
        update = 'Add a catastrophe';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = AddCatCoupled(state_x, state_y);
    case 14
        update='Delete a catastrophe';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = DelCatCoupled(state_x, state_y);
    case 15
        update = 'Resample catastrophes';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = ResampleCatastrophesCoupled(state_x, state_y);
    case 16
        update = 'Vary kappa';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = scaleCoupledKappa(state_x, state_y, mcmc);
    case 17
        error('Depracated');
        % update = 'Vary lambda';
        % [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        %     = scaleCoupledLambda(state_x, state_y, mcmc);
    case 18
        update='Move catastrophe to neighbour';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = MoveCatCoupled(state_x, state_y);
    case 19
        update = 'Vary XI for one leaf';
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = scaleCoupledXi(state_x, state_y, mcmc);
    case 21
        update = 'Vary beta';
        if ~VARYBETA && BORROWING
            disp('vary beta?')
            keyboard;
        end
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = scaleCoupledBeta(state_x, state_y, mcmc);
    otherwise
        error('Not implemented');
    end

    % Contribution to logq from catastrophe locations
    if OK_x && BORROWING && ismember(MV, [1:7, 11:12])
        logq_x = logq_x + catastropheScalingFactor(state_x, nstate_x);
    end
    if OK_y && BORROWING && ismember(MV, [1:7, 11:12])
        logq_y = logq_y + catastropheScalingFactor(state_y, nstate_y);
    end

    if OK_x
        [state_x, succ_x] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
            logq_x, state_x, nstate_x, OK_x, U_x, ignoreearlywarn);
    else
        succ_x = 0;
    end
    if OK_y
        [state_y, succ_y] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
            logq_y, state_y, nstate_y, OK_y, U_y, ignoreearlywarn);
    else
        succ_y = 0;
    end
end
