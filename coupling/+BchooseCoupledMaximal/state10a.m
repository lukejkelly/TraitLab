function [s, state] = state10a(cladeStatus)
    % Same leaves and clades as state10 but different trees
    % Used as a second tree in BchooseCoupledMaximalTest
    filePath = fullfile('coupling', '+BchooseCoupledMaximal', ...
                        sprintf('state10a-clades%s.mat', cladeStatus));
    state = getfield(load(filePath), 'state');
    s = state.tree;
end
