function [ x_T ] = solveOverTree( tEvents, tr, borPars )
% Function so solve a series of systems of DEs across a tree.
% Takes as inputs struct of tree events and transition rates.
% Returns the pattern means at the leaves.

% Pattern mean just before first branching event.
x_sminus = tr(1) / tr(2);

for k = 1:( size(tEvents, 2) - 1 )

    if tEvents(k).type == 1 % Branching event and solve.

        % Pass pattern means across branching event.
        x_s = transferMatrix( tEvents(k).L - 1, x_sminus, ...
            tEvents(k).loc, borPars );

    elseif tEvents(k).type == 2 % Catastrophe.

        % The duration is as per Ryder and Nicholls (2011).
        x_s = catastropheMeans( tEvents(k).L, x_sminus, ...
            -log( 1 - tr(4) ) / tr(2), tr, tEvents(k).loc, ...
            tEvents(k).K, borPars );

    elseif tEvents(k).type == 3 % Branch death.

        % Just assign pattern means and continue to solve as dead branches
        % are taken care of by patternMeans.
        x_s = x_sminus;

    end

    % Solve up to next event time.
    if (tEvents(k).time - tEvents(k + 1).time > 0)

        % Use pattern means approach for small L.
        x_tminus = patternMeans( tEvents(k).L, x_s, tEvents(k).time ...
            - tEvents(k + 1).time, tr, tEvents(k).K, borPars );

    else

        % In the unlikely event that two events occur simultaneously
        x_tminus = x_s;

    end

    % Set x_sminus to the current value and repeat loop.
    x_sminus = x_tminus;

end

% Return x_T.
x_T = x_sminus;

end
