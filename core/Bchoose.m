function [i, j, k, newage, logq, ncat, cat, loc] = Bchoose(state, mt, THETA, prior)

global ROOT WIDE OTHER MCMCCAT
N = 2 * state.NS - 1;
s = state.tree;

FAIL = 0;

i = ceil(N * rand);
while s(i).type == ROOT
   i = ceil(N * rand);
end
iP = s(i).parent;
iT = s(i).time;

if mt == WIDE

    if prior.isclade
        r=[];
        for a = 1:(2 * state.NS)
            if s(a).type < ROOT & length(s(a).clade) >= length(s(iP).clade) & length(s(s(a).parent).clade) >= length(s(iP).clade) & isempty(setdiff(s(iP).clade, s(a).clade)) & isempty(setdiff(s(s(a).parent).clade, s(iP).clade))
                r = [r, a];
            end
        end
        N = length(r);
    else
        r = 1:N;
    end

    if N>4
        j = r(ceil(N * rand));
        k = s(j).parent;
        while ( s(k).time <= iT || i == j || i == k )
            j = r(ceil(N * rand));
            k = s(j).parent;
        end
    else
        FAIL = 1;
        k = -1;
        j = -1;
    end

else

   if s(iP).type == ROOT
      FAIL = 1;
      k = -1;
      j = -1;
   else
      k = s(iP).parent;
      j = s(k).child(OTHER(s(iP).sibling));
   end

end

if FAIL || k == iP || j == iP
   newage = [];
   logq = -inf;
   [ncat, cat, loc, ~] = Bcats.failOutputs();
else
    [newage, logq_b] = Bchoose.sampleMarginal(i, j, k, s, THETA);
    if MCMCCAT
        [ncat, cat, loc, logq_c] = Bcats(state, i, j, k, newage);
    else
        [ncat, cat, loc, logq_c] = Bcats.failOutputs();
    end
    logq = logq_b + logq_c;
end
