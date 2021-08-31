function [newage_x, newage_y, logq_x, logq_y] ...
        = sampleCoupling(i, j_x, j_y, k_x, k_y, s_x, s_y, THETA)
    global ROOT

    % Get parameters
    if s_x(j_x).type == ROOT
        jT_x = s_x(j_x).time;
    else
        [new_minage_x, kT_x, logq_x] ...
            = BchooseCoupled.sampleCoupling.getUniformParameters(...
                i, j_x, k_x, s_x, THETA);
    end
    if s_y(j_y).type == ROOT
        jT_y = s_y(j_y).time;
    else
        [new_minage_y, kT_y, logq_y] ...
            = BchooseCoupled.sampleCoupling.getUniformParameters(...
                i, j_y, k_y, s_y, THETA);
    end

    % Sample new node ages
    if s_x(j_x).type == ROOT && s_y(j_y).type == ROOT
        % Both new ages exponential above respective current root
        [newage_x, newage_y] = maximalCouplingShiftedExponential(...
            jT_x, jT_y, THETA);

        logq_x = BchooseCoupled.sampleCoupling.getShiftedExponentialParameter(...
            i, jT_x, s_x, newage_x, THETA);
        logq_y = BchooseCoupled.sampleCoupling.getShiftedExponentialParameter(...
            i, jT_y, s_y, newage_y, THETA);
    elseif s_x(j_x).type == ROOT
        % Draw newage in x from shifted exponential and y from uniform
        [newage_y, newage_x] = maximalCouplingUniformShiftedExponential(...
            new_minage_y, kT_y, jT_x, THETA);

        logq_x = BchooseCoupled.sampleCoupling.getShiftedExponentialParameter(...
            i, jT_x, s_x, newage_x, THETA);
    elseif s_y(j_y).type == ROOT
        % Draw newage in x from uniform and y from shifted exponential
        [newage_x, newage_y] = maximalCouplingUniformShiftedExponential(...
            new_minage_x, kT_x, jT_y, THETA);

        logq_y = BchooseCoupled.sampleCoupling.getShiftedExponentialParameter(...
            i, jT_y, s_y, newage_y, THETA);
    else
        % Both newages uniform, logq depends on whether iP is the root or not
        [newage_x, newage_y] ...
            = maximalCouplingUniform(new_minage_x, kT_x, new_minage_y, kT_y);
    end
end
