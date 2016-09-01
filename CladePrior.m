function OK=CladePrior(prior,state,CHECKARG)

OK=1;
ncl=length(state.claderoot);
for k=1:ncl
    rn=state.claderoot(k);
    %rn is the lowest node covering all the nodes in clade k
    %state.tree(rn).unclade lists all the clades with the property that
    %node rn covers at least one leaf in the clade and at least one not in
    %the clade. If rn is the clade root for clade k then rn should not be
    %above any leaves which are not members of the clade.
    if prior.strongclades & ~isempty(state.tree(rn).unclade) & any(state.tree(rn).unclade==k)
        OK=0;
        if nargin==3
            disp(['clade ',prior.clade{k}.name,' is not a clade of the tree input to CladePrior()']);
            %keyboard; pause;
        else
            return;
        end
    end    
    rr=prior.clade{k}.rootrange;
    if ~isempty(rr)       
        if ~InRange(state.tree(rn).time,rr)
            OK=0;
            if nargin==3
                disp(['root of clade ',prior.clade{k}.name,' outside imposed range']);
                %keyboard; pause;
            else
                return;
            end
        end
    end
    ar=prior.clade{k}.adamrange;
    if ~isempty(ar)
        an=state.tree(state.claderoot(k)).parent;
        if ~InRange(state.tree(an).time,ar)
            OK=0;
            if nargin==3
                disp(['adam of clade ',prior.clade{k}.name,' outside imposed range']);
                %keyboard; pause;
            else
                return;
            end
        end
    end
end