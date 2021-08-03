function plotAll(d, n)
    clf;
    nexttile([3, 4])
    catsOnBranches.plotCDF(sum(d, 2), sum(n, 2));
    yyaxis left
    ylabel('(E)CDF');
    yyaxis right
    ylabel('Error');
    legend('ECDF', 'CDF', 'Error');
    title('Tree');
    for j = find(sum(d, 1))
        nexttile;
        catsOnBranches.plotCDF(d(:, j), n(:, j));
        title(sprintf('Branch %i', j));
    end
end
