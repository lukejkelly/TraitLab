function saveGlobal()

    % save('coupling/+MarkovCoupledTest/uncoupledVariables-20210603.mat', ...
    %      'mcmc', 'model', 'state_x', 'state_y');
    % save('coupling/+MarkovCoupledTest/coupledVariables-20210603.mat', ...
    %      'mcmc', 'model', 'state_x', 'state_y');


    globalVars = who('global');
    for iVar = 1:numel(globalVars)
        eval(sprintf('global %s', globalVars{iVar}));  % [EDITED]
    end

    save('coupling/+MarkovCoupledTest/globalVariables-20210603.mat', ...
         globalVars{:});
end
