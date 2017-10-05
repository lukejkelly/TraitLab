function[ x_t ] = transferMatrix(L, x_tminus, loc, borPars)
% A function which maps the pattern means, x_tminus, just before a
% branching event on branch loc (of L) to their corresponding values
% afterwards, x_t, when there are L + 1 branches.

% x_tminus is the pattern means before the split, loc the branching point
% (counting from the right-most branch in the tree representation in the
% plane) and L the number of languages.
% x_t is the pattern means after the split, i.e. at the beginning of the
% next period.

% Matrix of binary patterns before split and vector of their corresponding
% decimal representations.
P_de_tminus = ( 1:(2^L - 1) )';
P_bi_tminus = borPars(end).P_b(P_de_tminus, (end + 1 - L):end);

% Matrix of binary patterns after split and vector of their corresponding
% decimal representations.
P_bi_t = P_bi_tminus( :, [ 1:(L + 1 - loc), (L + 1 - loc):L ] );

% If a compiled version of the C function fastBi2De_c_par.c exists, we use
% it for the binary-to-decimal conversion, and, if not, the built-in Matlab
% function.
% P_de_t = fastBi2De_c_par(P_bi_t); % Binary-to-decimal function parallelised in C.
P_de_t = bi2de(P_bi_t, 'left-msb'); % Matlab's binary-to-decimal function.

% The P_de_t(i)th entry of x_t is the P_de_tminus(i)th entry of x_tminus.
% x_t is zero everywhere else.
x_t = zeros(2^(L + 1) - 1, 1);
x_t(P_de_t) = x_tminus;

end
