% Assume we're pasting in to Matlab running in the TraitLabSDLT-coupled folder
GlobalSwitches;
GlobalValues;

burnin = 5e2;
outPath = '../../Desktop/debug';
outFile_m = 'prior-m-0';
outFile_c = 'prior-c-0_x';

d_m = catsOnBranches.getBranchLengths(outPath, outFile_m, burnin);
d_c = catsOnBranches.getBranchLengths(outPath, outFile_c, burnin);

n_m = catsOnBranches.getCatastropheCounts(outPath, outFile_m, burnin);
n_c = catsOnBranches.getCatastropheCounts(outPath, outFile_c, burnin);

% Overall moments
[mObs_m, mExp_m, vObs_m, vExp_m] = catsOnBranches.getMoments(d_m, n_m);
[mObs_c, mExp_c, vObs_c, vExp_c] = catsOnBranches.getMoments(d_c, n_c);

disp([mObs_m, mExp_m, vObs_m, vExp_m]);
disp([mObs_c, mExp_c, vObs_c, vExp_c]);

% Plot distributions
FPropName = {'PaperPosition', 'PaperSize'};
[FW, FH] = deal(30, 25);
FPropVal = {[0, 0, FW, FH], [FW, FH]};

catsOnBranches.plotAll(d_m, n_m); pause(1);
set(gcf, FPropName, FPropVal); pause(1);
print(gcf, '-dpdf', sprintf('%s/%s.pdf', outPath, 'marginal')); pause(1);

catsOnBranches.plotAll(d_c, n_c); pause(1);
set(gcf, FPropName, FPropVal); pause(1);
print(gcf, '-dpdf', sprintf('%s/%s.pdf', outPath, 'coupled')); pause(1);
