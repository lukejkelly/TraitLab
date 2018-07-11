function [lPL] = logPredictiveFromStates(pathToData, testData, outFile, ...
                                         sInds, misDat, lostOnes)
% Function to compute the log predictive likelihood for a test data set using
% samples from the model fit to the corresponding test partition.
%   pathToData --- path to input nexus file directory
%   testData   --- name of test data nexus file (without .nex ending)
%   outFile    --- output file stem for experiments
%   sInds      --- indices of samples to use to compute log predictive
%   misDat     --- did we account for missing data in training data experiments
%   lostOnes   --- did we account for rare traits (discard singletons) in
%                  training data experiments.

  % Setting up
  GlobalSwitches; GlobalValues;
  global MISDAT LOSTONES;
  MISDAT = misDat; LOSTONES = lostOnes;

  % Clearing log-likelihood function
  clear('logLkd2_m');

  % Reading test data
  [~, cTest, ~, ~] = nexus2stype(sprintf('%s%s%s.nex', pathToData, ...
                                         filesep, testData));

  % Discarding any entries not satisfying the registration process
  if LOSTONES == 1
      rInds = find(sum(cTest.array == 1, 1) > 1);
  else
      rInds = find(any(cTest.array == 1, 1));
  end
  cTest.array = cTest.array(:, rInds);
  cTest.L = length(rInds);

  % Compute log predictive likelihood of test data
  lPL = zeros(length(sInds), 1);

  for k = 1:length(sInds)

    % Sampled parameters
    load(sprintf('saveStates%s%s-%05d.mat', filesep, outFile, sInds(k)), ...
         'state');

    % Replace training data in state by test data
    for i = 1:length(cTest.language)
      for j = state.leaves
        if strcmp(cTest.language(i), state.tree(j).Name)
          state.tree(j).dat = cTest.array(i, :);
        end
      end
    end
    state.L = cTest.L;

    % Compute log-likelihood
    clear('patternCounts');
    lPL(k) = logLkd2_m(state);

  end

end
