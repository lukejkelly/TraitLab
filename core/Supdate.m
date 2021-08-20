function [nstate, U, TOPOLOGY] = Supdate(state, i, newage)
    % mode node i to time newage
    % GKN; last modified by RJR on 09/05/07

    global LEAF ROOT ANST

    nstate = state;
    oldage = state.tree(i).time;
    nstate.tree(i).time = newage;

    switch state.tree(i).type
    case LEAF
        disp('Warning: Supdate should not be used to modify age of a leaf.');
        keyboard; pause;
        nstate.length = state.length + oldage - newage;
    case ANST
        nstate.length = state.length - oldage + newage;
    case ROOT
        nstate.length = state.length - 2 * oldage + 2 * newage;
    end

    U = above(i, state.tree, state.root);
    TOPOLOGY = 0;

end
