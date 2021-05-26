function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global BORROWING MCMCCAT;
    testCase.TestData.BORROWING = BORROWING;
    testCase.TestData.MCMCCAT = MCMCCAT;
    BORROWING = 0;
    MCMCCAT = 0;
    warning('Assume no borrowing or catastrophes');
end
