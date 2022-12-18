function dy = ode2BPkep( ~, y, mu, J2, R )
% ODE system for the two-body problem (Keplerian motion) in cartesian coordinates
%
% PROTOTYPE
% @(t,y) ode2BPkep(t, y, mu, 0, 0)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( a, e, i, bOmega, sOmega, theta ) [ L, ADIM, ANGLE ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% J2[1] J2 perturbation
% R[1] Radius of the planet, only relevant if J2 != 0
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L, ADIM, ANGLE ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
% Keplerian elements
a = y(1); e = y(2); i = y(3); bOmega = y(4); sOmega = y(5); theta = y(6);

% Intermediate
p = a*(1-e^2);
rnorm = p/(1+e*cos(theta));
h = sqrt(p*mu);

% Set the derivatives of the state
dy = [	0
		0
		0
		0
		0
		h/rnorm^2];
end