function [state,pa] = Markov(mcmc,model,state,ignoreearlywarn)
% [state,pa] = Markov(mcmc,model,state)
% written by GKN
% last modified by RJR on 07/05/07

global STOPRUN WIDE NARROW DEPNU VARYP ADAM VARYLAMBDA VARYMU TESTUP DONTMOVECATS BORROWING VARYBETA

if nargin<=3, ignoreearlywarn=0; end

acct=zeros(mcmc.update.Nmvs,1);
prop=zeros(mcmc.update.Nmvs,1);

for t=1:(mcmc.subsample)

    drawnow;
    STOPRUN = get(gcbo,'UserData');
    if STOPRUN
        %calculate proportion accepted for each update type
        pa=zeros(mcmc.update.Nmvs,1);
        pa(prop~=0)=acct(prop~=0)./prop(prop~=0);
        return
    end

    %choose an update type according to the cumulative dbn cmove
    r=rand;
    MV=find(r<mcmc.update.cmove,1,'first');
    prop(MV)=prop(MV)+1;

    %compute the candidate state
    OK=1;
    if MV==1
        update='RW node time between parent time and oldest child time';
        [i,newage,logq]=Schoose(state);
        [nstate,U,TOPOLOGY]=Supdate(state,i,newage);
    elseif MV==2
        update='Exchange nearest neighbours';
        [i,j,iP,jP,logq,OK]=Echoose(state,NARROW,model.prior);
        if OK
            [nstate,U,TOPOLOGY]=Eupdate(state,i,j,iP,jP);
        end
    elseif MV==3
        update='reconnect two edges randomly chosen across tree';
        [i,j,iP,jP,logq]=Echoose(state,WIDE,model.prior);
        if OK
            [nstate,U,TOPOLOGY]=Eupdate(state,i,j,iP,jP);
        end
    elseif MV==4
        update='Reconnect an edge into a nearby edge';
        [i,j,k,newage,logq]=Bchoose(state,NARROW,mcmc.update.theta,model.prior);
        OK=~isempty(newage);
        if OK
            [nstate,U,TOPOLOGY]=Bupdate(state,i,j,k,newage);
        end
    elseif MV==5
        update='Reconnect an edge into an edge chosen UAR over the tree';
        [i,j,k,newage,logq]=Bchoose(state,WIDE,mcmc.update.theta,model.prior);
        OK=~isempty(newage);
        if OK
            [nstate,U,TOPOLOGY]=Bupdate(state,i,j,k,newage);
        end
    elseif MV==6
        update='Rescale whole tree';
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        [nstate,U,TOPOLOGY,OK,logq]=Rscale(state,variation);
    elseif MV==7
        update='Rescale randomly chosen subtree';
        [nstate,U,TOPOLOGY,logq,OK]=RscaleSubTree(state,mcmc.update.del,mcmc.update.deldel);
    elseif MV==8
        update='Vary mu';
        if ~VARYMU, disp('vary mu ?'); keyboard;pause; end
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq=-log(variation);
        nstate=state;
        nstate.mu=state.mu*variation;
        TOPOLOGY=0;
        U=state.nodes;
        if DEPNU
            nstate.nu=state.nu/variation;
        end
        if DONTMOVECATS && nstate.mu<1e-5, OK=0; end
    elseif MV==9
        update='Vary p';
        if ~VARYP, disp('vary p ?'); keyboard;pause; end
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        pp=state.p^variation;
        logq=log(pp/(variation*state.p));
        nstate=state;
        nstate.p=pp;
        TOPOLOGY=0;
        U=state.nodes;
    elseif MV==10
        update='Missing data';
        if ~MISDAT, disp('reconstruct missing data ?'); keyboard;pause; end
        keyboard;
    elseif MV==11
        update='Vary leaves';
        vl=state.leaves(find([state.tree(state.leaves).leaf_has_timerange])); %#ok<FNDSB>
        OK=~isempty(vl);
        if OK
            nvl=size(vl,2);
            i=vl(ceil(rand*nvl));
            nstate=state;
            nstate.tree(i).time=state.tree(i).timerange(1)+rand*(state.tree(i).timerange(2)-state.tree(i).timerange(1));
            %nstate.tree(i).cat=(state.tree(i).cat-state.tree(i).time)/(state.tree(state.tree(i).parent).time-state.tree(i).time)*(nstate.tree(nstate.tree(i).parent).time -nstate.tree(i).time)+nstate.tree(i).time;
            %keyboard;
            %tr=reshape([state.tree(vl).timerange],2,size(vl,2));
            %cld=num2cell(tr(1,:)+diff(tr).*rand(1,size(vl,2)));
            %[nstate.tree(vl).time]=deal(cld{:});
            if nstate.tree(nstate.tree(i).parent).time>nstate.tree(i).time
                logq=0;
                TOPOLOGY=0;
                U=above(i,nstate.tree,nstate.root);
                nstate.length=state.length+state.tree(i).time-nstate.tree(i).time;
            else
                OK=0;
            end
        end
    elseif MV==12
        update='Rescale top tree';
        if ~model.prior.isclade
            error('Move 12 should not be selected if no clades');
        end
        [nstate,U,TOPOLOGY,logq,OK]=RscaleTopTree(state,model.prior,mcmc.update.del,mcmc.update.deldel);
    elseif MV==13
        update='Add a catastrophe';
        [nstate, U, OK, logq]=AddCat(state);
        TOPOLOGY=0;
    elseif MV==14
        update='Delete a catastrophe';
        [nstate, U, OK, logq]=DelCat(state);
        TOPOLOGY=0;
    elseif MV==15
        update='Vary rho';
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        nstate=state;
        nstate.rho=state.rho*variation;
        logq=-log(variation);
        TOPOLOGY=0;
        U=[];
    elseif MV==16
        update='Vary kappa';
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq=-log(variation); %GKN Feb 08 - was =0; was bug for kappa unif prior?
        nstate=state;
        if (variation*state.kappa <= 1) && (variation*state.kappa >= 0.25)
            nstate.kappa=variation*state.kappa;
        else
            OK=0;
        end
        TOPOLOGY=0;
        U=state.nodes;
        if OK && DEPNU
            nstate.nu=state.nu/variation;
        end
    elseif MV==17
        update='Vary lambda';
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq=-log(variation);
        nstate=state;
        nstate.lambda=state.lambda*variation;
        TOPOLOGY=0;
        U=state.nodes;
        if DEPNU
            nstate.nu=state.nu*variation;
        end
    elseif MV==18
        update='Move catastrophe to neighbour';
        [nstate,U,OK,logq]=MoveCat(state);
        TOPOLOGY=0;
    elseif MV==19
        %xi<1??
        update='Vary XI for one leaf';
        nstate=state;
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq=-log(variation);
        leaf=state.leaves(ceil(rand*length(state.leaves)));
        if variation*(1-state.tree(leaf).xi) <=1
            nstate.tree(leaf).xi=1-variation*(1-state.tree(leaf).xi);
        else
            OK=0;
        end
        U=above(leaf,state.tree,state.root);
        TOPOLOGY=0;
    elseif MV==20
        update='Vary XI for all leaves';
        nstate=state;
        V=[];
        variation=mcmc.update.del+rand*mcmc.update.deldel;
        logq = -2 * log(variation); % logq=2*log(variation);
        for i=state.leaves
            if state.tree(i).xi<1
                c=variation*(1-state.tree(i).xi);
                if c<=1
                    nstate.tree(i).xi=1-c;
                    logq = logq + log(variation); % logq=logq-log(variation);
                    V=[V,i];
                else
                    OK=0;
                    break;
                end
            end
        end
        U=above(V,state.tree,state.root);
        TOPOLOGY=0;
    elseif MV == 21 % LUKE
        update = 'Vary beta';
        % if ~VARYBETA && BORROWING, disp('vary beta?'), keyboard; pause; end
        variation = mcmc.update.del + rand * mcmc.update.deldel;
        logq = -log(variation);
        nstate = state;
        nstate.beta = state.beta * variation;
        TOPOLOGY = 0;
        U = state.nodes;
    end % End of move section.

    % Contribution to logq from catastrophe locations
    if OK && BORROWING && ismember(MV, [1:7, 11:12])
        logq = logq + catastropheScalingFactor(state, nstate);
    end

    if OK && model.prior.isclade
        if TOPOLOGY
            nstate=UpdateClades(nstate,U,size(model.prior.clade,2));
        end
        OK=CladePrior(model.prior,nstate);
    end

    if OK

        if ~all(sort(find([state.tree.type] == 0)) == sort(find([nstate.tree.type] == 0)))
            disp('Leaf variables mismatch')
            keyboard; return;
        end

%         if TOPOLOGY
%             nstate.length=TreeLength(nstate.tree,nstate.root);
%         end

        %disp(update);
        if TOPOLOGY && DONTMOVECATS
            nstate.cat(:)=1;
            nstate.cat(nstate.root)=0;
            nstate.cat(nstate.tree(nstate.root).parent)=0;
            nstate.ncat=sum(nstate.cat);
        end
        [nstate, OK]=MarkRcurs(nstate,U,TOPOLOGY,ignoreearlywarn);

        % Likelihood calculations.
        if BORROWING
            [nstate.loglkd, nstate.fullloglkd] = logLkd2(nstate); % LUKE
        else
            nstate.loglkd     = LogLkd(nstate);
            nstate.fullloglkd = LogLkd(nstate, nstate.lambda);
        end

        % Log-prior for tree and catastrophes.
        nstate.logprior = LogPrior(model.prior, nstate);

        % Log-prior for parameters: mu, beta, etc.
        logpp = LogPriorParm(state, nstate);

        % Log-Hastings ratio in acceptance step.
        logr = nstate.logprior - state.logprior + nstate.loglkd - state.loglkd + logq + logpp;

        % Acceptance step.
        %GKN Jan 08 to handle big tree, llkd is -inf so do random walk (3rd OR)
        if ( (logr>0) || (log(rand)<logr) || isinf(state.loglkd) )
            state=nstate;
            acct(MV)=acct(MV)+1;
        end
        if TESTUP && check(state,[])
            disp(['Error in update:',update]);
            check(nstate,[])
            keyboard;pause;
        end
    end

end

% % % % % %update lambda|g,mu - assumes 1/lambda prior
% % % % % if ~VARYLAMBDA
% % % % %     state.lambda=Lambda(state);
% % % % %     if DEPNU
% % % % %         state.nu=state.kappa*state.lambda/state.mu;
% % % % %     end
% % % % % end

%calculate proportion accepted for each update type
%pa=zeros(mcmc.update.Nmvs,1);
%pa(prop~=0)=acct(prop~=0)./prop(prop~=0);
pa=acct./prop; %Changed to display NaN when the update was not proposed. RJR 19/03/09
