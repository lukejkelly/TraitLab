% Fast replacements for Communications Toolbox functions de2bi and bi2de.
% Run from main TraitLab folder after adding borrowing/ to the path.
% Include -fopenmp flags if using parallelised versions.

% [Parallelised] Binary-to-decimal function.
mex -outdir borrowing -output bi2de...
    borrowing/fastBi2De.c ...
    CFLAGS="\$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic" ... % -fopenmp
    LDFLAGS="\$LDFLAGS" % -fopenmp

% [Parallelised] Decimal-to-binary function.
mex -outdir borrowing -output de2bi ...
    borrowing/fastDe2Bi.c ...
    CFLAGS="\$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic" ... % -fopenmp
    LDFLAGS="\$LDFLAGS" % -fopenmp
