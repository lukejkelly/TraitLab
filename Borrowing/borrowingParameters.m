function [ borPars ] = borrowingParameters( L )
% Function to return the parameters necessary for indexing patterns and the
% non-catastrophe differential system in the borrowing model on a tree with
% L leaves. For catastrophes and branch deaths, we compute the parameters
% as required.

% Creating an empty struct to hold the parameters.
borPars = struct( 'P_b', cell(L, 1), 'S', cell(L, 1), 'D', cell(L, 1), ...
    'A_d', cell(L, 1), 'A_b', cell(L, 1), 'b', cell(L, 1) );

% Indices for populating arrays.
inds = (1:(2^L - 1))';

% Populating P_b array. We only need to compute P_b for L as P_b for k < L
% is just a sub-array.
% If a compiled version of the C function fastDe2Bi_c_par exists, we use it
% and if not, the slower Matlab function.
% borPars(L).P_b = fastDe2Bi_c_par(inds, L);
borPars(L).P_b = de2bi(inds, L, 'left-msb');

% Populating S. Once again, we only need S for L.
borPars(L).S = sum(borPars(L).P_b, 2);

% Populating D. likewise, we only need D for L.
% borPars(L).D = repmat(inds, 1, L) ...
%     + (-1).^borPars(L).P_b * diag(2.^((L - 1):(-1):0));
borPars(L).D = bsxfun(@plus, inds, ...
    (1 - 2 * borPars(L).P_b) * diag(2.^((L - 1):(-1):0)));

% The entries in the rate matrix A depend on the death and borrowing
% parameters. In which case, we shall compute two separate matrices, A_d
% and A_b, and pass them to the patternMeans function which can scale them
% by the corresponding parameters before adding them together.
for l = 2:L
    
    % Number and index of patterns.
    eta = 2^l - 1;
    inds = 1:eta;

    % Pattern sizes.
    S = borPars(L).S(inds);
    
    % Generating the indices of patterns which differ from the elements of
    % P_b by one entry.
    D = borPars(L).D( inds, (L + 1 - l):end );

    % Setting up differential system.
    
    % Vectors to populate A_d and A_b.
    dInds = find(D);
    iInds = 1 + mod(dInds - 1, eta);
    jInds = D(dInds);

    A_i  = [(1:eta)'         ; iInds];
    A_j  = [(1:eta)'         ; jInds];
    A_sd = [-S               ;                 (S(jInds) > S(iInds))];
    A_sb = [-S .* (1 - S / l); S(jInds) / l .* (S(jInds) < S(iInds))];

    % Constructing sparse rate matrices. Matrices are stored in CSC (compressed
    % sparse column) format, so for our purposes we store the transpose.
    borPars(l).A_d = sparse( A_j, A_i, A_sd );
    borPars(l).A_b = sparse( A_j, A_i, A_sb );

    % Vector of birth rates.
    borPars(l).b = sparse( S(inds) == 1 );
    
end

% Converting to lower-precision integer types in order to save memory and
% speed up computation. We can't work with smaller types to start with as
% the MTIMES (*) operation in the creaton of D is not supported for
% integer types.
if L < 32; borPars(L).D = uint32( borPars(L).D ); end

end
