function Y = timed2BP( r, v, mu, ngrid, time, nPeriods )
% Simple orbit determination, for fast plotting (low accuracy)
%
% PROTOTYPE
% Y = timed2BP( r, v, mu, ngrid )
%
% INPUT:
% r[3x1] Position ( rx, ry, rz ) [ L ]
% v[3x1] Velocity ( vx, vy, vz ) [ L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% ngrid[1] Number of points (low number is recommended)
% time[1] Length of time plotted. If not specified, it reverts to 1 orbital period
%
% OUTPUT:
% Y[ngridx3] Position of ngrid points in the specified time / orbital period [ L ]
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

% Solver options
options = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-7 );

if (nargin == 5)
	if isscalar(time)
		T = linspace( 0, time, ngrid )';
	else
		T = linspace( 0, time(1), ceil(ngrid/2) )';
		[ ~, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu, 0, 0), T, y0, options );
		T = linspace( 0, time(2), ceil(ngrid/2) )';
		[ ~, Y2 ] = ode113( @(t,y) ode_2bp(t,y,mu, 0, 0), T, y0, options );
		Y = [flip(Y1);Y2];
		return
	end
else
	% Calculate semi-major axis and orbit period
	rNorm = vecnorm(r);
	vNorm = vecnorm(v);
	a = mu/(2*mu/rNorm-vNorm^2);
	Torb = 2*pi*sqrt( a^3/mu );

	if (nargin == 4)
		% 1 period time grid
		T = linspace( 0, Torb, ngrid )';
	else 
		T = linspace( 0, nPeriods*Torb, ngrid )';
	end
end

% Integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu, 0, 0), T, y0, options );
Y = Y(:,1:3);
end