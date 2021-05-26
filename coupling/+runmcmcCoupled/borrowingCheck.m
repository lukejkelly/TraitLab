% Borrowing check inside MCMC loop
function borrowingCheck(state)
    % Luke 10/02/14
    % Checking to make sure catastrophe locs aren't getting screwed around.
    % We compare the number of catastrophe locations with the number we expect.
    global BORROWING
    if BORROWING
        for i = 1:(2 * state.NS)
            if state.cat(i) ~= length(state.tree(i).catloc)
                sprintf('Catastrophe mismatch on <pa(%d), %d>', i, i)
                keyboard; pause;
            end
        end
    end
end
