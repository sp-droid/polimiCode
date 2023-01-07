function [Y,T] = timed2BP( y0, mu, opts, ngrid, time, nPeriods )
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
%		• show[bool]		-  Print options at the end
%		• perturbShow[bool] -  Print perturbation trigger states
%		• TinPeriods[bool]	-  Divide time grid by orbital period
%		• AbsTol[1]			-  Solver absolute tolerance, defaults to 1e-7
%		• RelTol[1]			-  Solver relative tolerance, defaults to 1e-6
%		• solver[@func]		-  Solver used, defaults to @ode113
%		• keplerian[bool]	-  If true, use keplerian elements (only J2 perturbation is implemented with kep. elems.)
%		• Rearth[1]			-  Earth radius
%		• muSun[1]		    -  Sun's gravitational parameter
%		• muMoon[1]		    -  Moon's gravitational parameter
%		• sunPos[@func]		-  rSun = @(t) func(t, args)
%		• moonPos[@func]	-  rMoon = @(t) func(t, args)
%		• J2[1]				-  Perturbation J2 term
% ngrid[1]		- Number of points (low number is recommended)
%
% PERTURBATIONS TRIGGERS:
% Simple J2     - Rearth, J2
% Sun           - muSun, sunPos
% Moon          - muMoon, moonPos
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
% Y[ngridx1]	- Time grid [ T ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
%% Options
if not(isstruct(opts))
	opts = struct;
end
if ~isfield(opts,'show')
	opts.show = false;
end
if ~isfield(opts,'perturbShow')
	opts.perturbShow = false;
end
if ~isfield(opts,'TinPeriods')
	opts.TinPeriods = false;
end

%% Perturbations
perturbs = struct;

% J2 perturbation
if isfield(opts,'j2Pert')
	perturbs.j2Pert = true;
else
	perturbs.j2Pert = false;
end
% 3rd body perturbation: Sun
if isfield(opts,'sunThirdBody')
	perturbs.sunThirdBody = true;
else
	perturbs.sunThirdBody = false;
end
% 3rd body perturbation: Moon
if isfield(opts,'moonThirdBody')
	perturbs.moonThirdBody = true;
else
	perturbs.moonThirdBody = false;
end
% Earth gravitational model 1996
if isfield(opts,'egm96') && isfield(opts,'wEarth')
	perturbs.egm96 = true;
else
	perturbs.egm96 = false;
end
% Relativistic effects
if isfield(opts,'relativEffect')
	perturbs.relativEffect = true;
else
	perturbs.relativEffect = false;
end
% Atmospheric drag
if isfield(opts,'drag')
	perturbs.drag = true;
else
	perturbs.drag = false;
end
% Solar radiation pressure
if isfield(opts,'srp')
	perturbs.srp = true;
else
	perturbs.srp = false;
end

%% Define spatial integration parameters
% Keplerian or cartesian 
if ~isfield(opts,'keplerian')
	opts.keplerian = false;
end
if opts.keplerian
	a = y0(1);
else
	rNorm = vecnorm(y0(1:3));
	vNorm = vecnorm(y0(4:6));
	a = mu/(2*mu/rNorm-vNorm^2);
end

%% Define time integration parameters
% Calculate orbit period
Torb = 2*pi*sqrt( a^3/mu );

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

% Define t0
if (interval(2)>interval(1))
    opts.t0 = interval(1);
else
    opts.t0 = interval(2);
end

%% Solver
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

if opts.keplerian
	odeSystem = @(t,y) ode2BPkep(t, y, mu, opts, perturbs);
else
	odeSystem = @(t,y) ode2BPcar(t, y, mu, opts, perturbs);
end

%% Integration
if (interval(1)==0)
	T = linspace( 0, interval(2), ngrid )';
	[ ~, Y ] = opts.solver( odeSystem, T, y0, solverOpts );
	if (interval(2)<0)
		T = flip(T);
		Y = flip(Y);
	end
else
	T1 = linspace( 0, interval(1), ceil(ngrid/2) )';
	[ ~, Y1 ] = opts.solver( odeSystem, T1, y0, solverOpts );
	T2 = linspace( 0, interval(2), ceil(ngrid/2) )';
	[ ~, Y2 ] = opts.solver( odeSystem, T2, y0, solverOpts );
	T = [flip(T1);T2];
	Y = [flip(Y1);Y2];
end

%% Wrap angles in case we are working with keplerian elements
if opts.keplerian
	Y(:,3) = wrapTo2Pi(2*Y(:,3))/2;	% i		 [0, pi]
	Y(:,4) = wrapTo2Pi(Y(:,4)); 	% bOmega [0, 2pi]
	Y(:,5) = wrapTo2Pi(Y(:,5)); 	% sOmega [0, 2pi]
	Y(:,6) = wrapTo2Pi(Y(:,6)); 	% theta  [0, 2pi]
else
%% Change time grid unit to orbital periods, if requested
if opts.TinPeriods
	T = T/Torb;
end

%% Show options
if opts.show
	opts = orderfields(opts);
	display(opts)
end
%% Show perturbation trigger states
if opts.perturbShow
	display(perturbs)
end
end