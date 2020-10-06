function [ obs, miss, rl_c ] = observedPatternCounts( state )
% Function return the pattern counts for the leaf order given by
% [state.tree(i).ord == L, ..., state.tree(i).ord == 1]. We shall
% manipulate this vector of counts to calculate the likelihood (instead of
% manipulating the leaf order of the sampled tree and pattern means, etc.).
% The variables returned are n_To and n_Tm representing the number of
% observations of each fully and partially observed pattern.

% Global variables.
global LEAF

% Number of leaves.
L = size(state.tree, 2) / 2; % state.NS;

% Tree.
s = state.tree;

% Choosing an arbitrary right-left order based on leaf labels. The vector
% rl_c contains the leaf-indices with the first entry corresponding to the
% index of the right-most leaf, etc.
rl_c = find([s.type] == LEAF);

% Initialising data array: (# leaves) x (# patterns).
P_b = zeros( L, size(s( rl_c(1) ).dat, 2) );

% Populating data array.
for k = 1:L

    % First row corresponds to right-most leaf.
    P_b(k, :) = s( rl_c(k) ).dat;

end

% Discarding 'zero' patterns.
P_b(:, ~any(P_b)) = [];

% Separating fully and partially-observed entries.
P_bo = P_b(:, ~any(P_b == 2) );
P_bm = P_b(:, any(P_b == 2) );

% Fully-observed patterns.
if ~isempty(P_bo)

    % The entries in the count vector correspond to patterns 00···01, ...,
    % 11···11 when the right-left order of the leaves is given by rl.

    % Reorientating P_bo.
    P_bo = fliplr( P_bo' );

    % Decimal representation of fully observed patterns (if there are any).
    % We use parallelised C function if its compiled version exists.
	% P_do = fastBi2De_c_par(P_bo);
    P_do = bi2de(P_bo, 'left-msb');

    % Initialising vector of counts.
    n_To = zeros(2^L - 1, 1);

    % Populating the count vector.
    for i = P_do'

        n_To(i) = n_To(i) + 1;

    end

    % Creating a struct with the binary representations of observed
    % patterns and their respective counts.
    inds = find(n_To);

    % For the decimal-to-binary conversion, we use a C function if its
    % compiled version is available.
	% obs.pattern = fastDe2Bi_c_par(inds, L);
	obs.pattern = de2bi(inds, L, 'left-msb');
    obs.count = n_To(inds);

else

    % If there are no fully-observed patterns then we set obs to 0. We
    % can't set it to [] as this will pose a problem when we use persistent
    % variables.
    obs = 0;

end

% Partially-observed patterns.
if ~isempty(P_bm)

    % For each partially-observed pattern, we require the pattern it refers
    % to, the corresponding vector of missing indices and matrix of
    % possible true underlying patterns.

    % Reorienting P_bm.
    P_bm = fliplr( P_bm' );

    % Unique entries - we only need the counts of the patterns we
    % (partially) observed.
    P_bmu = unique(P_bm, 'rows', 'stable');

    % Creating a struct which, for each partially observed pattern contains
    % the pattern, the count and a list of possible true underlying
    % patterns and a vector indicating.
    miss = struct('pattern', cell(size(P_bmu, 1), 1), ...
                  'count', cell(size(P_bmu, 1), 1), ...
                  'underlying', cell(size(P_bmu, 1), 1));

    % Populating the count vector.
    for i = 1:size(P_bmu, 1)

        % Recording the pattern with missing data.
        miss(i).pattern = P_bmu(i, :);

        % Number of observations of partially-observed pattern P_bmu(i, :).
        miss(i).count = sum( all( bsxfun( @eq, P_bmu(i, :), P_bm ), 2 ) );

        % Indices of entries with data missing.
        m_inds = find( P_bmu(i, :) == 2 );

        % Generating possible underlying patterns.
        % Using C function, if its compiled version exists.
        % perms = fastDe2Bi_c_par( 0:(2^length(m_inds) - 1), length(m_inds) );
        % Using Matlab function otherwise.
        perms = de2bi( 0:(2^length(m_inds) - 1), length(m_inds), 'left-msb' );

        p_poss = repmat( P_bmu(i, :), length(perms), 1 );
        p_poss(:, m_inds) = perms;

        % removing any patterns consisting solely of zeros and creating the
        % corresponding struct entry.
        miss(i).underlying = p_poss(any(p_poss, 2), :);

    end

else

    % If there are no partially-observed patterns then we set miss to 0. We
    % can't set it to [] as this will pose a problem when we use persistent
    % variables.
    miss = 0;

end

end
