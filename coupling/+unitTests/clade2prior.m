function prior = clade2prior(clade, rootmax)
    % Derived from fullsetup>initModel
    global FLAT YULE;
    prior.clade = clade;
    prior.strongclades = 1;
    n = size(prior.clade, 2);
    prior.isupboundadamclade = logical([]);
    prior.upboundclade=[];
    prior.isclade = ~isempty(clade);
    if rand < 0.5
        prior.type = FLAT;
    else
        prior.type = YULE;
    end
    prior.rootmax = rootmax;
    for k = 1:n
        if (~isempty(prior.clade{k}.adamrange) && ~isempty(prior.clade{k}.rootrange))
            if (prior.clade{k}.adamrange(1) < prior.clade{k}.rootrange(1))
                prior.clade{k}.adamrange(1) = prior.clade{k}.rootrange(1);
                warning(['Setting the lower bound on origin in clade ', prior.clade{k}.name, ' equal to the higher lower bound on root.']);
            end
            if (prior.clade{k}.rootrange(2) > prior.clade{k}.adamrange(2))
                prior.clade{k}.rootrange(2) = prior.clade{k}.adamrange(2);
                warning(['Setting the upper bound on root in clade ', prior.clade{k}.name, ' equal to the lower upper bound on origin.']);
            end
        end
        if (~isempty(prior.clade{k}.adamrange) && ((prior.type == FLAT && prior.clade{k}.adamrange(2) < prior.rootmax) || (prior.type == YULE && ~isinf(prior.clade{k}.adamrange(2)))))
            prior.upboundclade=[prior.upboundclade, k];
            prior.isupboundadamclade=[prior.isupboundadamclade, true];
            prior.clade{k}.lowlim = prior.clade{k}.adamrange(1);
        else
            if (~isempty(prior.clade{k}.rootrange) && ((prior.type == FLAT && prior.clade{k}.rootrange(2) < prior.rootmax)) || (prior.type == YULE && ~isinf(prior.clade{k}.rootrange(2))))
                prior.upboundclade=[prior.upboundclade, k];
                prior.isupboundadamclade=[prior.isupboundadamclade, false];
                prior.clade{k}.lowlim = prior.clade{k}.rootrange(1);
            end
        end
    end
end
