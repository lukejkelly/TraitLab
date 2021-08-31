function [i_x, i_y, j_x, j_y, iP_x, iP_y, jP_x, jP_y, logq_x, logq_y, OK_x, ...
          OK_y] = EchooseCoupled(state_x, state_y, mt, prior)
    % Sample from a maximal coupling of uniform distributions over valid (i, j)
    % pairs in states x and y
    global NARROW WIDE

    switch mt
    case NARROW
        v_x = EchooseCoupled.narrowTargets(state_x);
        v_y = EchooseCoupled.narrowTargets(state_y);

        [i_x, i_y] = discreteUniformCoupledSample(v_x, v_y);

        [j_x, iP_x, jP_x, OK_x] = EchooseCoupled.narrowOutputs(i_x, state_x, prior);
        [j_y, iP_y, jP_y, OK_y] = EchooseCoupled.narrowOutputs(i_y, state_y, prior);
    case WIDE
        v_x = EchooseCoupled.wideTargets(state_x, prior);
        v_y = EchooseCoupled.wideTargets(state_y, prior);

        if ~isempty(v_x) && ~isempty(v_y)
            [h_x, h_y] = discreteUniformCoupledSample(v_x, v_y);
        elseif ~isempty(v_x)
            h_x = discreteUniformSample(v_x);
            h_y = [];
        elseif ~isempty(v_y)
            h_x = [];
            h_y = discreteUniformSample(v_y);
        else
            h_x = [];
            h_y = [];
        end
        [i_x, j_x, iP_x, jP_x, OK_x] = EchooseCoupled.wideOutputs(state_x, h_x);
        [i_y, j_y, iP_y, jP_y, OK_y] = EchooseCoupled.wideOutputs(state_y, h_y);
    end

    logq_x = 0;
    logq_y = 0;
end
