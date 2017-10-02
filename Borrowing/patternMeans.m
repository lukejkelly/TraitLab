function [x_tminus] = patternMeans(L, x_s, T, tr, K, borPars)
% Routine to compute pattern means / expected pattern frequencies across an
% interval of length T between branching points on a tree.
% There are currently L branches.
% x_s are the pattern means at the start of the interval and x_tminus are
% those at the end.
% K are the indices of dead branches as these must be accounted for.
% tr = [lambda, mu, beta, kappa] transition rates and catastrophe severity
% (not required for this function).

% Creating persistent variable to hold parameters for differential equation
% from the previous iteration when there are offset leaves. If the indices
% of offset leaves don't change from one state to the next, we reuse those
% from the previous step. If the indices do change, we permute the previous
% versions.
persistent offsetPars

if isempty(offsetPars)

    % Number of leaves on tree.
    M = size(borPars, 1);

    % Creating an empty struct to hold parameters.
    offsetPars = struct( 'A_d', cell(M), 'A_b', cell(M), 'b', cell(M), ...
        'K', cell(M) );

end

% Number of patterns across L branches.
eta = 2^L - 1;
inds = 1:eta;

% Setting up differential system.
if isempty(K)

    % Creating A and b matrices for input to differential equation solver.
    A = borPars(L).A_d * tr(2) + borPars(L).A_b * tr(3);
    b = borPars(L).b * tr(1);

else

    % Number of dead branches.
    ndead = size(K, 2);

    % If offsetPars(L, ndead).K is empty then we need to create
    % corresponding A and b parameters.
    if isempty( offsetPars(L, ndead).K )

        % Update offsetPars(L, ndead).K.
        offsetPars(L, ndead).K = K;

        % Pattern sets.
        P_b = borPars(end).P_b(inds, (end + 1 - L):end);

        % Pattern sizes.
        S = borPars(end).S(inds);

        % Indices of patterns which differ from rows of P_b by one entry.
        D = borPars(end).D(inds, (end + 1 - L):end);

        % Indices of alive branches.
        inds_a    = 1:L;
        inds_a(K) = [];

        % Number of evolving entries in each pattern.
        S_a = sum( P_b(:, L + 1 - inds_a), 2 );

        % Removing from D the entries corresponding to dead branches.
        D(:, L + 1 - K) = 0;

        % Number of extant branches.
        L_a = L - ndead;

        % Setting up differential system.

        % Empty vectors with which we shall populate A_d and A_b.

% This section has been replaced by the somewhat harder to read code below.
%         % Vectors of (i, j) indices.
%         A_i = zeros( eta * (L_a + 1) - L_a, 1 );
%         A_j = zeros( eta * (L_a + 1) - L_a, 1 );
%
%         % Vectors of A(i, j) entries.
%         A_sd = zeros( eta * (L_a + 1) - L_a, 1 );
%         A_sb = zeros( eta * (L_a + 1) - L_a, 1 );
%
%         % Populating A_i, A_j, A_sd and A_sb.
%         k = 1;
%
%         for i = inds
%
%             % Diagonal entries before scaling by \mu and \beta.
%             % A(i, i) = A_d(i, i) + A_b(i, i),
%             %         = -S_a(i) * (1 + (L_a - S_a(i))/L_a.
%             A_i(k) = i;
%             A_j(k) = i;
%             A_sd(k) = - S_a(i);
%             A_sb(k) = - S_a(i) * ( L_a - S_a(i) ) / L_a;
%
%             k = k + 1;
%
%             % Indices of off-diagonal entries corresponding to entries in
%             % the row vector A(i, :). A(i, j_inds) is the sub-vector of
%             % indices of patterns which flow into pattern i.
%             j_inds = D( i, D(i, :) > 0 );
%             len = length( j_inds );
%             k_inds = k:(k + len - 1);
%
%             % Populating off-diagonal entries.
%             A_i( k_inds ) = i;
%             A_j( k_inds ) = j_inds;
%
%             % A_sd(i, j) = 1 if j \in S_p^+.
%             A_sd( k_inds ) = ( S_a(j_inds) > S_a(i) );
%
%             % A_sd(i, j) = S_a / L_a if j in S_p^-.
%             A_sb( k_inds ) = S_a(j_inds) / L_a .* ( S_a(j_inds) < S_a(i) );
%
%             % Updating index parameter.
%             k = k + len;
%
%         end

        dInds = find(D);
        iInds = 1 + mod(dInds - 1, eta);
        jInds = D(dInds);

        A_i  = [(1:eta)'               ; iInds];
        A_j  = [(1:eta)'               ; double(jInds)];
        A_sd = [-S_a                   ;                     (S_a(jInds) > S_a(iInds))];
        A_sb = [-S_a .* (1 - S_a / L_a); S_a(jInds) / L_a .* (S_a(jInds) < S_a(iInds))];

        % Constructing sparse matrices (we store the transpose).
        offsetPars(L, ndead).A_d = sparse(A_j, A_i, A_sd);
        offsetPars(L, ndead).A_b = sparse(A_j, A_i, A_sb);

        % Vector b of birth rates.
        offsetPars(L, ndead).b = sparse( (S == 1 & S_a == 1) );

    % If offsetPars(L, ndead).K is non-empty, is it the same as the current
    % K? If so, then we can reuse it. If not, then we permute the previous
    % A and b matrices to obtain those for the current set of indices of
    % dead branches.
    elseif ~all( sort( offsetPars(L, ndead).K ) == sort(K) )

        % Sorted current and old indices of dead branches.
        K_old  = sort( offsetPars(L, ndead).K );
        K_curr = sort(K);

        % Patterns to be permuted.
        P_b = borPars(end).P_b(inds, (end + 1 - L):end);

        % Efficiently finding indices on which K_old and K_curr differ.

        % Empty vectors.
        v_old  = zeros(L, 1);
        v_curr = zeros(L, 1);

        % Populating vectors.
        v_old(K_old)   = -1;
        v_curr(K_curr) = 1;

        % Entries with a -1 need to be switched with those with a +1.
        v_upd  = v_old + v_curr;
        i_old  = find(v_upd == -1);
        i_curr = find(v_upd == 1);

        % Permutation operation.
        P_b(:, L + 1 - [i_old, i_curr]) = P_b(:, L + 1 - [i_curr, i_old]);

        % Decimal representations of permuted patterns.
        % If compiled parallelised C function exists, we use it for binary-to-
        % decimal operation. And if not, the built-in Matlab function.
        % inds_curr = fastBi2De_c_par(P_b);
        inds_curr = bi2de(P_b, 'left-msb');

        % Permuting entries in offsetPars(L).A_d/A_b/b to create the
        % matrices we require for the ODE.
        A_d = offsetPars(L, ndead).A_d(inds_curr, inds_curr);
        A_b = offsetPars(L, ndead).A_b(inds_curr, inds_curr);
        b   = offsetPars(L, ndead).b(inds_curr);

        % Updating offsetPars(L, ndead) fields A_d, A_b and b to be the
        % permuted matrices we now require for K_curr.
        offsetPars(L, ndead).A_d = A_d;
        offsetPars(L, ndead).A_b = A_b;
        offsetPars(L, ndead).b   = b;

        % Updating offsetPars(L, ndead) field K to reflect current offset
        % lineages.
        offsetPars(L, ndead).K = K_curr;

    end

    % Forming parameters of differential equation.
    A = offsetPars(L, ndead).A_d * tr(2) + offsetPars(L, ndead).A_b * tr(3);
    b = offsetPars(L, ndead).b * tr(1);

end

% Solution of ODE.
x = ode45( @(t, x) A' * x + b, [0, T], x_s );

% Output
x_tminus = x.y(:, end);

end
