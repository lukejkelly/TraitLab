% Fast replacements for Communications Toolbox functions de2bi and bi2de.
% Run from main TraitLab folder after adding Borrowing/ to the path.
% Include -fopenmp flags if using parallelised versions.

% [Parallelised] Binary-to-decimal function.
mex -IBorrowing -outdir Borrowing -output bi2de...
    Borrowing/fastBi2De.c ...
    CFLAGS="\$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic" ... % -fopenmp
    LDFLAGS="\$LDFLAGS" % -fopenmp

% [Parallelised] Decimal-to-binary function.
mex -IBorrowing -outdir Borrowing -output de2bi ...
    Borrowing/fastDe2Bi.c ...
    CFLAGS="\$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic" ... % -fopenmp
    LDFLAGS="\$LDFLAGS" % -fopenmp
