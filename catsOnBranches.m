% Assume we're pasting in to Matlab running in the TraitLabSDLT-coupled folder
GlobalSwitches;
GlobalValues;

burnin = 1e3;
outPath = '../../Desktop/debug';
outFile_m = 'sim-n-m-0';
outFile_c = 'sim-n-c-0_x';

d_m = catsOnBranches.getBranchLengths(outPath, outFile_m);
d_c = catsOnBranches.getBranchLengths(outPath, outFile_c);

n_m = catsOnBranches.getCatastropheCounts(outPath, outFile_m, burnin);
n_c = catsOnBranches.getCatastropheCounts(outPath, outFile_c, burnin);

% Overall moments
[mObs_m, mExp_m, vObs_m, vExp_m] = catsOnBranches.getMoments(d_m, n_m);
[mObs_c, mExp_c, vObs_c, vExp_c] = catsOnBranches.getMoments(d_c, n_c);

disp([mObs_m, mExp_m, vObs_m, vExp_m]);
disp([mObs_c, mExp_c, vObs_c, vExp_c]);

% Plot distributions
catsOnBranches.plotAll(d_m, n_m);
exportgraphics(gcf, sprintf('%s/%s.pdf', outPath, 'marginal'));

catsOnBranches.plotAll(d_c, n_c);
exportgraphics(gcf, sprintf('%s/%s.pdf', outPath, 'coupled'));
