function res = LogPriorParm(state, nstate)

% calculates the log of the ratio between the priors for the parameters mu,
% lambda, nu, kappa, rho, beta in state and nstate.

global VARYNU VARYBETA

% mu and beta priors are Gamma(a, a)
a = 1e-3;

res = lpdfG(nstate.mu, a, a) - lpdfG(state.mu, a, a);

% LJK 06/21, if VARYRHO/VARYLAMBDA then integrated out analytically
% if VARYRHO, res = res + LogRhoPrior(nstate.rho) - LogRhoPrior(state.rho); end
% if VARYLAMBDA, res = res + log(state.lambda) - log(nstate.lambda); end

if VARYNU
    res = res + log(state.nu) - log(nstate.nu);
end

if VARYBETA
    res = res + lpdfG(nstate.beta, a, a) - lpdfG(state.beta, a, a);
end

%if nstate.kappa<0.1, res=-inf; end

%if nstate.kappa>.9, res=-inf; end %This line because of machine precision
