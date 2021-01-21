function compareDistributions(testCase, xObs, xExp)
    [~, edges] = histcounts([xObs, xExp], 50);
    nObs = histcounts(xObs, edges, 'Normalization', 'cdf');
    nExp = histcounts(xExp, edges, 'Normalization', 'cdf');
    plot(edges', [zeros(2, 1), [nObs; nExp]]', ':', 'LineWidth', 2);
    legend('Obs', 'Exp', 'Location', 'SouthEast');
    axis('tight');
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end
