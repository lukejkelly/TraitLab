function [rf, lf, lb] = getDistributions(catc, ncat, wc, wn)
    % Handles to simulate and evaluate log-probability for forward move, and
    % log-probability for return move
    rf = @() mnrnd(ncat, wn);
    lf = @(z) log(mnpdf(z(:)', wn));
    lb = log(mnpdf(catc(:)', wc));
end
