function[ x_tminus ] = catastropheMeans(L, x_s, T, tr, loc, K, borPars)
% Function to compute the pattern means as a result of a catastrophe on
% lineage loc (counting from right to left).
% Takes the same inputs as patternMeans with the addition of loc:
%  -- L: number of lineages
%  -- x_s: vector of initial expected pattern frequencies
%  -- T: length of interval
%  -- tr: rate parameters
%  -- loc: catastrophe branch
%  -- K: offset leaves on interval, if any
%  -- borPars: parameters of ODEs, etc.
% Outputs pattern means at the end of the catastrophe period, which is
% a period of virtual time during which deaths on and borrowings in to
% loc occur.

% We find pairs of patterns which differ by one at location loc then form
% the corresponding rate matrix. We generate the pairs of patterns instead
% of finding them.

% K lists the indices of branches which have died. Dead branches don't
% affect which patterns get paired up during a catastrophe (due to the fact
% that paired patterns only differ at the index corresponding to the
% catastrophe branch. We alter the borrowing rates between patterns during
% a catastrophe by the number of extant branches and evolving traits in the
% source pattern.

% Updated by LJK on 20/11/2014 to use a for-loop and the explicit solution
% of 2 x 2  ODEs rather than one large system of ODEs and ODE solver.

% Number of patterns.
inds = 1:(2^L - 1);

% Generating the indices of patterns which only differ by their entry on
% the catastrophe branch.
D = borPars(end).D(inds, end + 1 - loc);

% Set of possible patterns across L branches.
P_b = borPars(end).P_b(inds, (end + 1 - L):end);

% Parameters which depend on whether there are offset leaves or not.
if isempty(K)

    % There are no dead branches; that is, K is empty.

    % Number of evolving traits in each pattern.
    S_a = borPars(end).S(inds);

    % Number of extant branches.
    L_a = L;

else

    % List of dead branches is non-empty.

    % Indices of extant branches.
    inds_a = 1:L;
    inds_a(K) = [];

    % Number of evolving traits in each pattern.
    S_a = sum( P_b(:, L + 1 - inds_a), 2 );

    % Number of extant branches.
    L_a = L - length(K);

end

% Updating pattern means.

% Empty vector of pattern means.
x_tminus = zeros(2^L - 1, 1);

% There is one pattern with a 1 on loc and 0s elsewhere which arise due to
% new traits being born and existing traits generating the pattern dying
% out.

% Index of pattern...
inds_b = 2^(loc - 1);

% ... and updating its mean.
x_tminus(inds_b) = exp(-tr(2) * T) * x_s(inds_b) ...
    + (1 - exp(-tr(2) * T)) * tr(1) / tr(2);

% Updating the rest of the patterns. For each pattern with a 0 on the
% catastrophe branch, we update it along with its pair which is identical
% on the off-catastrophe-branch entries.

% Indices of patterns which have a 0 on the catastrophe branch.
inds_c = inds( ~P_b(:, L + 1 - loc) );

% We have that \dot(x)(t) = A * x(t), where x(t) = [x_pt, x_qt], s(p) < s(q),
% m = \mu, b = \beta · s(p) / L_a and A = [-b, m; b, -m]. It is thus the case
% that x(t) = expm(A * t) x(0), where expm(A * T) ) =
%   / m + b exp(-T(b + m))      m - m exp(-T(b + m)) \
%   | ----------------------, ---------------------- |
%   |          b + m                   b + m         |
%   |                                                |.
%   | b - b exp(-T(b + m))      b + m exp(-T(b + m)) |
%   | ----------------------, ---------------------- |
%   \          b + m                   b + m         /
% We precompute a1 and a2 where
%   a1 = (m + b exp(-T(b + m))) / (b + m),
%   a2 = (m - m exp(-T(b + m))) / (b + m),
% for various b to avoid repeated exponentiation in update functions.

% Setting mu, b for each possible pattern size...
mu = tr(2);
b = (0:L_a)' * tr(3) / L_a;

% ... and computing a1 and a2.
a1 = (mu + b  .* exp(-T * (b + mu))) ./ (b + mu);
a2 = (mu - mu .* exp(-T * (b + mu))) ./ (b + mu);

% Updating the pattern means.
x_tminus(inds_c) = a1(S_a(inds_c) + 1) .* x_s(inds_c)  ... % Contribution from x_p0.
    + a2(S_a(inds_c) + 1) .* x_s(D(inds_c)); % Contribution from x_q0.

% Mass is conserved so we find x_qt as the difference of the initial values
% and x_pt.
x_tminus(D(inds_c)) = x_s(inds_c) + x_s(D(inds_c)) - x_tminus(inds_c);

end
