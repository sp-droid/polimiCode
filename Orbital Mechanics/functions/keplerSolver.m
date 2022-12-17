function E = keplerSolver(t, e, a, mu, t0, E0)
%
% PROTOTYPE
% E = keplerSolver(t, e, a, mu, t0, E0)
%
% INPUT:
% t[1] Time vector
% e[1] Eccentricity
% a[1] Semi-major axis
% mu[1] Gravitational parameter
% t0[1] Reference time
% E0[1] Reference eccentric anomaly
%
% OUTPUT:
% E[1] Eccentric anomaly
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-09-21: v1
%
% -------------------------------------------------------------------------
% Mean motion and initial guess for E
n = sqrt(mu/a^3);
Eguess = n*t+e*sin(n*t)/(1-sin(n*t+e)+sin(n*t));

% Define function to solve for
func = @(E) E-e*sin(E)-n*t;

% Solve
E = fzero(func, Eguess);
end