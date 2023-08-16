function [state, succ] = MarkovCoupledCommon(...
        mcmc, model, state, ignoreearlywarn, MV, u_mh)
    global BORROWING

    OK = 1;
    switch MV
    case 11
        update='Vary leaves';
        vl=state.leaves(find([state.tree(state.leaves).leaf_has_timerange]));
        OK=~isempty(vl);
        if OK
            nvl=size(vl,2);
            i=vl(ceil(rand*nvl));
            nstate=state;
            if isinf(state.tree(i).timerange(2))
                % LK 08/23 if clade time is ignored then upper bound is inf so
                % sample uniformly between 0 and root time instead
                nstate.tree(i).time = state.tree(i).timerange(1) + rand * (state.tree(state.root).time - state.tree(i).timerange(1));
            else
                nstate.tree(i).time = state.tree(i).timerange(1) + rand * (state.tree(i).timerange(2) - state.tree(i).timerange(1));
            end
            if nstate.tree(nstate.tree(i).parent).time>nstate.tree(i).time
                logq=0;
                TOPOLOGY=0;
                U=above(i,nstate.tree,nstate.root);
                nstate.length=state.length+state.tree(i).time-nstate.tree(i).time;
            else
                OK=0;
            end
        end
    case 20
        update='Vary XI for all leaves';
        nstate=state;
        V=[];
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq = -2 * log(variation);
        for i=state.leaves
            if state.tree(i).xi<1
                c=variation*(1-state.tree(i).xi);
                if c<=1
                    nstate.tree(i).xi=1-c;
                    logq = logq + log(variation);
                    V=[V,i];
                else
                    OK=0;
                    break;
                end
            end
        end
        U=above(V,state.tree,state.root);
        TOPOLOGY=0;
    otherwise
        error('Move %d is depracated or does not exist', MV);
    end

    % Contribution to logq from catastrophe locations
    if OK && BORROWING && ismember(MV, [1:7, 11:12])
        logq = logq + catastropheScalingFactor(state, nstate);
    end

    if OK
        [state, succ] = MarkovUpdateState(update, model, u_mh, TOPOLOGY, ...
            logq, state, nstate, OK, U, ignoreearlywarn);
    else
        succ = 0;
    end
end
