function res=LogPriorParm(state,nstate)

% res=LogPriorParm(state,nstate) calculates the log of the ratio between
% the priors for the parameters mu, lambda, nu, kappa and rho in state and
% nstate.

global VARYNU VARYLAMBDA VARYBETA % BORROWING

%if BORROWING
    % We can't use a prior of 1/mu as the posterior is unbounded so we use
    % a Gamma(0.001, 0.001) prior with high variance to approximate it.
    res = -0.999 * (log(nstate.mu) - log(state.mu)) - 0.001 * (nstate.mu - state.mu);
%else
%    res=log(state.mu)-log(nstate.mu);
%end

% LJK 06/21, if VARYRHO then it's integrated out analytically, fixed otherwise
% if VARYRHO, res=res+LogRhoPrior(nstate.rho)-LogRhoPrior(state.rho);end
if VARYLAMBDA, res=res+log(state.lambda)-log(nstate.lambda); end
if VARYNU, res=res+log(state.nu)-log(nstate.nu); end

if VARYBETA
    % As with mu, we use a Gamma(0.001, 0.001) prior for beta in place of
    % 1/beta to ensure the posterior is proper.
    res = res + (-0.999 * (log(nstate.beta) - log(state.beta)) ...
                 - 0.001 * (nstate.beta - state.beta));
end

%if nstate.kappa<0.1, res=-inf; end

%if nstate.kappa>.9, res=-inf; end %This line because of machine precision
