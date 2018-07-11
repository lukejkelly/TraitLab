function [] = partitionData(path, data, pTrain)
% Randomly partition path/file.nex into separate training and test sets with
% pTrain of the data points going into the training partition.
%   path   --- directory containing nexus file we want to split up
%   data   --- name of the nexus file (without the .nex ending)
%   pTrain --- percentage of data (on average) to go into training partition,
%              set to 2/3 if undefined.

  % Setting up
  GlobalSwitches; GlobalValues;
  [sFull, contentFull, ~, cladeFull] ...
    = nexus2stype(sprintf('%s%s%s.nex', path, filesep, data));

  % Split into pTrain training and 1 - pTrain test
  if isempty(pTrain); pTrain = 2 / 3; end
  iTrain = (rand(1, contentFull.L) < pTrain);
  iTest  = ~iTrain;

  [sTrain, sTest] = deal(sFull);
  for l = find([sFull.type] == LEAF)
    sTrain(l).dat(iTest) = [];
    sTest(l).dat(iTrain) = [];
  end

  % Write new data files
  fidTrain = fopen(sprintf('%s%s%s-train.nex', path, filesep, data), 'w+');
  fprintf(fidTrain, stype2nexus(sTrain, 'SIM_B train portion', 'BOTH', [], ...
                                cladeFull));
  fclose(fidTrain);

  fidTest = fopen(sprintf('%s%s%s-test.nex', path, filesep, data), 'w+');
  fprintf(fidTest, stype2nexus(sTest, 'SIM_B test portion', 'BOTH', [], ...
                                cladeFull));
  fclose(fidTest);

end
