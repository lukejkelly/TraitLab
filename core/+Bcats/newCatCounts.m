function [ncat, cat] = newCatCounts(state, i, j, jn)
    % New counts on nodes h and p
    root = state.root;
    [p, ~, h] = Bcats.pqh(state, i);

    pc = state.cat(p);
    hc = state.cat(h);
    jc = state.cat(j);

    if p == root
        % h is new root so loses catastrophes
        hn = 0;
        pn = jc - jn;
        ncat = state.ncat - hc;
    elseif j == root
        % p is new root, j gains jn catastrophes
        hn = hc + pc;
        pn = 0;
        ncat = state.ncat + jn;
    else
        % total number of catastrophes remains constant
        hn = hc + pc;
        pn = jc - jn;
        ncat = state.ncat;
    end
    cat = struct('j', jn, 'h', hn, 'p', pn);
end
