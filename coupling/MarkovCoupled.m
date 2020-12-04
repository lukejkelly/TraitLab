function [state_x, state_y, pa_x, pa_y] = MarkovCoupled(mcmc, model, ...
    state_x, state_y, ignoreearlywarn)

    global STOPRUN

    if nargin<=4, ignoreearlywarn=0; end

    [acct_x, acct_y] = deal(zeros(mcmc.update.Nmvs, 1));
    prop=zeros(mcmc.update.Nmvs,1);

    for t=1:(mcmc.subsample)

        drawnow;
        STOPRUN = get(gcbo,'UserData');
        if STOPRUN
            pa_x = acct_x ./ prop;
            pa_y = acct_y ./ prop;
            return
        end

        % Common roles across trees have common indices
        state_y = housekeeping(state_x, state_y);
        state_y.tree = superHousekeeping(state_x.tree, state_y.tree);
        % For now we assume that state_x.tree was only ever a valid modification
        % of state_y so leaves and Adam have the same indices in both, as do all
        % nodes within a clade

        % MCMC acceptance probability
        u_mh = rand;

        %choose an update type according to the cumulative dbn cmove
        r=rand;
        MV=find(r<mcmc.update.cmove,1,'first');
        prop(MV)=prop(MV)+1;

        % MCMC updates
        if ismember(MV, [1, 4, 5, 8, 15, 16, 17, 19, 21])
            % Attempt maximal coupling
            [state_x, succ_x, state_y, succ_y] = MarkovCoupledMaximal(...
                mcmc, model, state_x, state_y, ignoreearlywarn, MV, u_mh);
            % if succ_x && succ_y
            %     fprintf('Successful coupling on %i\n', MV);
            % end
        else
            % Common random number coupling
            rng_state = rng;
            [state_x, succ_x] = MarkovCoupledCommon(mcmc, model, state_x, ...
                ignoreearlywarn, MV, u_mh);
            rng(rng_state);
            [state_y, succ_y] = MarkovCoupledCommon(mcmc, model, state_y, ...
                ignoreearlywarn, MV, u_mh);
            % Uncoupled chains may make different calls to rng
            rng('shuffle');
        end
        acct_x(MV) = acct_x(MV) + succ_x;
        acct_y(MV) = acct_y(MV) + succ_y;
    end
    pa_x = acct_x ./ prop;
    pa_y = acct_y ./ prop;
end
