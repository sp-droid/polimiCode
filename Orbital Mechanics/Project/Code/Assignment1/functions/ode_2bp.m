function dy = ode_2bp (~, y, mu)
% ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
% 
% VERSIONS:
%   2022-09-20
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

% Position and velocity
r = y(1:3);
v = y(4:6);

rnorm = norm(r);


% Derivative od state

dy = [ v; (-mu/rnorm^3)*r];
end




