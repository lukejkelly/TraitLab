function[ D ] = simBorCatDeath(tEvents, tr)
% Function to simulate patterns on a catastrophe tree with borrowing and
% branch deaths.
% -- tEvents is the output of stype2Events detailing branching, death and
%      catastrophe events, etc.
% -- tr is the vector of transition rates and catastrophe severity (lambda, mu
%      beta, kappa).
% During a catastrophe, the usual process is simulated but only events on
% the catastrophe branch are accepted.

% Parameters
lambda = tr(1);
mu     = tr(2);
beta   = tr(3);
kappa  = tr(4);

% Patterns just after branching event at root.
D = ones( poissrnd(lambda / mu), 1);

for k = 1:( size(tEvents, 2) - 1 )

    % Total number of branches
  	L = tEvents(k).L;

    % Total number of branches still alive.
    L_a = L - length( tEvents(k).K );

    % Indices of alive and dead branches.
    a_inds = 1:L;
    a_inds( L + 1 - tEvents(k).K ) = [];
    d_inds = L + 1 - tEvents(k).K;

    if tEvents(k).type == 1 % Branching event.

        % Branching location.
        loc = tEvents(k).loc;

        % Pass D across branching event.
        D = D( :, [ 1:( L - loc), (L - loc):(L - 1) ] );

    elseif tEvents(k).type == 2 % Catastrophe event.

        % Some parameters.
        len = - log(1 - kappa) / mu;
        cat_loc = L + 1 - tEvents(k).loc; % Counting from right.
        t = 0;

        while t < len

            % Time to next event is exponential. We only simulate births
            % and deaths on and transfers in to the catastrophe branch.
            % Transfers are restricted to come from branches which are
            % still extant.
            tlambda = lambda; % Birth rate on cat_loc
            tmu = mu * sum( D(:, cat_loc) ); % Likewise deaths.
            tbeta = beta * sum( sum( D(:, a_inds) ) ); % Borrow from extant branches.

            dt = exprnd( 1 / (tlambda + tmu + tbeta) );

            if t + dt < len

                % Event type.
                p = [tlambda, tmu, tbeta] / (tlambda + tmu + tbeta);

                etype = find( mnrnd(1, p) );

                if etype == 1 % Pattern birth.

                    % Add on correct pattern of size 1 to D.
                    v = zeros(1, L);
                    v(cat_loc) = 1;

                    D = [ D; v ];

                elseif etype == 2 % Trait death.

                    % Find indices of traits that are present.
                    inds = find( D(:, cat_loc) == 1 );

                    % Pick one uniformly at random and kill it.
                    D( inds( ceil(size(inds, 1) * rand) ), cat_loc ) = 0;

                else % Borrowing event.

                    % Select a pattern with a weight proportional to the
                    % number of alive patterns.
                    r_ind = mnrnd( 1, sum(D(:, a_inds), 2) / ...
                        sum( sum( D(:, a_inds) ) ) );

                    % Select a slot to transfer to.
                    c_ind = mnrnd(1, ones(1, L_a) / L_a);

                    % If selected slot is cat_loc, carry out move.
                    if a_inds( logical(c_ind) ) == cat_loc

                        D( logical(r_ind), cat_loc ) = 1;

                    end

                end % End if.

            end % End if.

            % Advance time.
            t = t + dt;

        end % End while.

    end % End if.

    % Some parameters.
    len = tEvents(k).time - tEvents(k + 1).time;
    t = 0;

    while t < len

        % Time to next event is exponential.
        tlambda = lambda * L_a;
        tmu = mu * sum( sum( D(:, a_inds) ) );
        tbeta = beta * sum( sum( D(:, a_inds) ) );

        dt = exprnd( 1 / (tlambda + tmu + tbeta) );

        if t + dt < len

            % Event type.
            p = [tlambda, tmu, tbeta] / (tlambda + tmu + tbeta);

            etype = find( mnrnd(1, p) );

            if etype == 1 % Pattern birth.

                % Choose a branch and create the pattern of size one.
                v = zeros(1, L);
                v( a_inds( logical( mnrnd(1, ones(1, L_a) / L_a) ) ) ) = 1;

                % Add on a randomly chosen pattern of size 1 to D.
                D = [ D; v ];

            elseif etype == 2 % Trait death.

                % Select a pattern with a weight proportional to the
                % number of alive patterns.
                r_ind = logical( mnrnd( 1, sum(D(:, a_inds), 2) / ...
                    sum( sum( D(:, a_inds) ) ) ) );

                % Select an entry to kill.
                c_ind = logical( mnrnd(1, D(r_ind, a_inds) / ...
                    sum( D(r_ind, a_inds) ) ) );

                % Set the corresponding entry in D to zero.
                D(r_ind, a_inds(c_ind) ) = 0;

            else % Borrowing event.

                % Select a pattern with a weight proportional to the
                % number of alive patterns to be the source of the
                % borrowing event.
                r_ind = logical( mnrnd( 1, sum(D(:, a_inds), 2) / ...
                    sum( sum( D(:, a_inds) ) ) ) );

                % Select a slot to transfer to.
                c_ind = logical( mnrnd(1, ones(1, L_a) / L_a) );

                % Set the corresponding entry in D to 1.
                D(r_ind, a_inds(c_ind) ) = 1;

            end % End if.

        end % End if.

        % Advance time.
        t = t + dt;

    end % End while.

end
