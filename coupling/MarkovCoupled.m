function [state_x, state_y, pa_x, pa_y, pa_xy] = MarkovCoupled(mcmc, model, ...
    state_x, state_y, ignoreearlywarn)

    global STOPRUN

    if nargin<=4, ignoreearlywarn=0; end

    [acct_x, acct_y, acct_xy] = deal(zeros(mcmc.update.Nmvs, 1));
    prop=zeros(mcmc.update.Nmvs,1);

    coupledBefore = checkCoupling(state_x, state_y);

    for t=1:(mcmc.subsample)

        drawnow;
        STOPRUN = get(gcbo,'UserData');
        if STOPRUN
            pa_x = acct_x ./ prop;
            pa_y = acct_y ./ prop;
            return
        end

        % MCMC acceptance probability
        u_mh = rand;

        %choose an update type according to the cumulative dbn cmove
        r=rand;
        MV=find(r<mcmc.update.cmove,1,'first');
        prop(MV)=prop(MV)+1;

        % MCMC updates
        if ismember(MV, [1, 4, 5, 6, 7, 8, 15, 16, 17, 19, 21])
            % Attempt maximal coupling
            [state_x, succ_x, state_y, succ_y] = MarkovCoupledMaximal(...
                mcmc, model, state_x, state_y, ignoreearlywarn, MV, u_mh);
        else
            % Common random number coupling
            rng_state = rng;
            [state_x, succ_x] = MarkovCoupledCommon(mcmc, model, state_x, ...
                ignoreearlywarn, MV, u_mh);
            rng(rng_state);
            [state_y, succ_y] = MarkovCoupledCommon(mcmc, model, state_y, ...
                ignoreearlywarn, MV, u_mh);
            % Uncoupled chains may make different number of calls to rng
            rng('shuffle');
        end

        % TODO: remove move uncoupling check after further testing
        coupledAfter = checkCoupling(state_x, state_y);
        if coupledBefore && ~coupledAfter
            save(sprintf('uncoupling move %s.mat', datetime));
            error('uncoupling at %s due to move', datetime);
        else
            coupledBefore = coupledAfter;
        end
        state_y = housekeeping(state_x, state_y);
        % TODO: remove housekeeping coupling check after further testing
        if coupledBefore && ~checkCoupling(state_x, state_y)
            save(sprintf('uncoupling housekeeping %s.mat', datetime));
            error('uncoupling at %s due to housekeeping', datetime);
        end
        acct_x(MV) = acct_x(MV) + succ_x;
        acct_y(MV) = acct_y(MV) + succ_y;
        acct_xy(MV) = acct_xy(MV) + succ_x * succ_y;
    end
    pa_x = acct_x ./ prop;
    pa_y = acct_y ./ prop;
    pa_xy = acct_xy ./ prop;
end
