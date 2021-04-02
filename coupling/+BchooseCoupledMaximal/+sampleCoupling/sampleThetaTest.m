function tests = sampleThetaTest
    % Compare to what Bchoose does
    tests = functiontests(localfunctions);
end

function checkDistributionTest(testCase)
    thetaObs = nan(1e4, 1);
    for i = 1:length(thetaObs)
        thetaObs(i) = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
    end
    log10ThetaExp = -3 + rand(size(thetaObs)) * 2;

    histogram(log10(thetaObs), linspace(-3, -1, 51), 'Normalization', 'pdf');
    hold on;
    histogram(log10ThetaExp,  linspace(-3, -1, 51), 'Normalization', 'pdf');
    hold off;
    refline(0, 1 / 2);
    legend('Obs', 'Exp', 'Location', 'SouthEast');
    axis('tight');

    v = input('Do the histograms in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end
