function Y = onePeriod2BP( r, v, mu, ngrid )
% Simple orbit determination, for fast plotting (low accuracy)
%
% PROTOTYPE
% Y = onePeriod2BP( r, v, mu, ngrid )
%
% INPUT:
% r[3x1] Position ( rx, ry, rz ) [ L ]
% v[3x1] Velocity ( vx, vy, vz ) [ L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% ngrid[1] Number of points (low number is recommended)
%
% OUTPUT:
% Y[ngridx3] Position of ngrid points for 1 orbital period [ L ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% State vector
y0 = [r; v];

% Calculate semi-major axis and orbit period
rNorm = vecnorm(r);
vNorm = vecnorm(v);
a = mu/(2*mu/rNorm-vNorm^2);
Torb = 2*pi*sqrt( a^3/mu );

% 1 period time grid
T = linspace( 0, Torb, ngrid )';

% Solver options
options = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-7 );

% Integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu, 0, 0), T, y0, options );
Y = Y(:,1:3);
end