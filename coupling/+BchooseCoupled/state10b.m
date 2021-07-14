function [s, state] = state10b(cladeStatus)
    % Same leaves and clades as state10a but different trees
    % Used as a second tree in BchooseCoupledTest
    filePath = fullfile('coupling', '+BchooseCoupled', ...
                        sprintf('state10b-clades%s.mat', cladeStatus));
    state = getfield(load(filePath), 'state');
    s = state.tree;
end
