function dy = ode_2bp( ~, y, mu, J2, R )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% J2[1] J2 perturbation
% R[1] Radius of the planet, only relevant if J2 != 0
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);

% aJ2 term
aJ2 = 5*r(3)^2/rnorm^2;
aJ2 = [r(1)/rnorm*(aJ2-1)
    r(2)/rnorm*(aJ2-1)
    r(3)/rnorm*(aJ2-3)];
aJ2 = aJ2*1.5*J2*mu*R^2/rnorm^4;

% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r+aJ2 ];
end