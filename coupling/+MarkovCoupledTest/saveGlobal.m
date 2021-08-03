function saveGlobal()

    % save('coupling/+MarkovCoupledTest/uncoupledVariables-20210803.mat', ...
    %      'mcmc', 'model', 'state_x', 'state_y');
    % state_y = state_x;
    % save('coupling/+MarkovCoupledTest/coupledVariables-20210803.mat', ...
    %      'mcmc', 'model', 'state_x', 'state_y');

    % save('coupling/+MarkovCoupledTest/uncoupledPriorVariables-20210803.mat', ...
    %      'mcmc', 'model', 'state_x', 'state_y');

    globalVars = who('global');
    for iVar = 1:numel(globalVars)
        eval(sprintf('global %s', globalVars{iVar}));
    end

    save('coupling/+MarkovCoupledTest/globalPriorVariables-20210803.mat', ...
         globalVars{:});
end
