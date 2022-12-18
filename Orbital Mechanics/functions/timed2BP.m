function Y = timed2BP( y0, mu, opts, ngrid, time, nPeriods )
% Orbit forward and/or backward propagation
%
% SIMPLEST PROTOTYPE
% Y = timed2BP( [r0;v0], mu, [], ngrid )
%
% INPUTS:
% y0			- State vector:
%		• y0[6carx1]			- Cartesian ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
%		• y0[6kepx1]			- Keplerian ( a, e, i, bOmega, sOmega, theta) [ L, ADIM, ANGLE]. Necessary to set opts.keplerian to True
% mu[1]			- Gravitational parameter of the primary [L^3/T^2]
% opts[struc]	- Structure with several options:
%		• AbsTol[1]			-  Solver absolute tolerance, defaults to 1e-7
%		• keplerian[bool]	-  True if using keplerian elements, defaults to false or cartesian
%		• RelTol[1]			-  Solver relative tolerance, defaults to 1e-6
%		• show[bool]		-  Print options at the end, defaults to false
%		• solver[@func]		-  Solver used, defaults to @ode113
% ngrid[1]		- Number of points (low number is recommended)
%
% OPTIONAL INPUTS:
% time			- Time interval in seconds:
%		• omitted			- Defaults to linspace [0, 1 period]
%		• time[1 (>0)]		- Translates into linspace [0, time]
%		• time[1 (<0)]		- Translates into linspace [time, 0]
%		• time[2x1]			- Translates into linspace [time(1), time(2)]
% nPeriods		- Time interval in orbital periods, which if defined superseeds time:
%		• nPeriods[1 (>0)]	- Translates into linspace [0, nPeriods]
%		• nPeriods[1 (<0)]	- Translates into linspace [nPeriods, 0]
%		• nPeriods[2x1]		- Translates into linspace [nPeriods(1), nPeriods(2)]
%
% OUTPUTS:
% Y[ngridx3]	- Position of ngrid points in the specified time / orbital period [ L ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------
%% Options
if not(isstruct(opts))
	opts = struct;
end
if ~isfield(opts,'show')
	opts.show = false;
end

%% Define spatial integration parameters
% Keplerian or cartesian 
if ~isfield(opts,'keplerian')
	opts.keplerian = false;
end
if opts.keplerian
	odeSystem = @(t,y) ode2BPkep(t, y, mu, 0, 0);
else
	odeSystem = @(t,y) ode2BPcar(t, y, mu, 0, 0);
end

% Solver used
if ~isfield(opts,'solver')
	opts.solver = @ode113;
end

% Solver relative tolerance
solverOpts = odeset();
if ~isfield(opts,'RelTol')
	opts.RelTol = 1e-6;
end
solverOpts.RelTol = opts.RelTol;

% Solver absolute tolerance
if ~isfield(opts,'AbsTol')
	opts.AbsTol = 1e-7;
end
solverOpts.AbsTol = opts.AbsTol;

%% Define time integration parameters
% Using time
if (nargin == 5)
	% [0, time] or [time, 0]
	if isscalar(time)
		interval = [0;time];
	% [time(1), time(2)]
	else
		interval = time;
	end
% Using nPeriods
else						
	% Calculate semi-major axis and orbit period
	rNorm = vecnorm(y0(1:3));
	vNorm = vecnorm(y0(4:6));
	a = mu/(2*mu/rNorm-vNorm^2);
	Torb = 2*pi*sqrt( a^3/mu );

	% [0, 1 period]
	if (nargin == 4)
		% 1 period time grid
		interval = [0;Torb];
	else
		% [0, nPeriods] or [nPeriods, 0]
		if isscalar(nPeriods)
			interval = [0;nPeriods*Torb];
		% [nPeriods(1), nPeriods(2)]
		else
			interval = Torb*nPeriods;
		end
		
	end
end

%% Integration
if (interval(1)==0)
	T = linspace( 0, interval(2), ngrid )';
	[ ~, Y ] = opts.solver( odeSystem, T, y0, solverOpts );
else
	T = linspace( 0, interval(1), ceil(ngrid/2) )';
	[ ~, Y1 ] = opts.solver( odeSystem, T, y0, solverOpts );
	T = linspace( 0, interval(2), ceil(ngrid/2) )';
	[ ~, Y2 ] = opts.solver( odeSystem, T, y0, solverOpts );
	Y = [flip(Y1);Y2];
end

%% Show options
if opts.show
	opts = orderfields(opts);
	disp(opts)
end
end