% Fast replacements for Communications Toolbox functions de2bi and bi2de.
% Run from main TraitLab folder after adding borrowing/ to the path.
% Include -fopenmp flags if using parallelised versions.

% Depending on your platform, you may need to replace single quotes with double 
% https://fr.mathworks.com/help/matlab/ref/mex.html#mw_663100fe-c5a5-4a18-9171-7ab49fe7d509
% One may also need to replace $ by \$

% [Parallelised] Binary-to-decimal function.
mex -outdir borrowing -output bi2de...
    borrowing/fastBi2De.c ...
    CFLAGS='$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic' ... % -fopenmp
    LDFLAGS='$LDFLAGS' % -fopenmp

% [Parallelised] Decimal-to-binary function.
mex -outdir borrowing -output de2bi ...
    borrowing/fastDe2Bi.c ...
    CFLAGS='$CFLAGS -O3 -std=c11 -march=native -Wall -pedantic' ... % -fopenmp
    LDFLAGS='$LDFLAGS' % -fopenmp
