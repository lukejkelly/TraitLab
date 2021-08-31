function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global BORROWING MCMCCAT;
    testCase.TestData.BORROWING = BORROWING;
    testCase.TestData.MCMCCAT = MCMCCAT;
end
