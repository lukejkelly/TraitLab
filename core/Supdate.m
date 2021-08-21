function [nstate, U, TOPOLOGY] = Supdate(state, i, newage, cat, loc)
    % mode node i to time newage, update catastrophes
    % GKN; last modified by RJR on 09/05/07, LJK on 21/08/21

    global LEAF ROOT ANST MCMCCAT

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

    if MCMCCAT
        nstate = Supdate.moveCats(nstate, i, cat, loc);
        U = above(state.tree(i).child, state.tree, state.root);
    else
        U = above(i, state.tree, state.root);
    end
    TOPOLOGY = 0;

end
