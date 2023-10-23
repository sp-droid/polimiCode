function dy = ode2BPcar( t, y, mu, opts, perturbs )
% ODE system for the two-body problem (Keplerian motion) in cartesian coordinates
%
% PROTOTYPE
% @(t,y) ode2BPcar(t, y, mu, 0, 0)
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
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3); v = y(4:6);

% Intermediate
rnorm = norm(r);

% Primary gravity acceleration
aGrav = (-mu/rnorm^3)*r;

aXYZ = [0;0;0];
if perturbs.j2Pert
	aXYZ = aXYZ + opts.j2Pert(r);
end
if perturbs.sunThirdBody
    aXYZ = aXYZ + opts.sunThirdBody(r,t);
end
if perturbs.moonThirdBody
    aXYZ = aXYZ + opts.moonThirdBody(r,t);
end
if perturbs.egm96
    % t-t0 assures the initial time of the propagation is 0
    % opts.theta0 is necessary to get the initial rotation of the Earth
    thetaG = opts.wEarth*(t-opts.t0)+opts.theta0;
    aXYZ = aXYZ + opts.egm96(r, thetaG);
end
if perturbs.relativEffect
    aXYZ = aXYZ + opts.relativEffect(r,v);
end
if perturbs.drag
    aXYZ = aXYZ + opts.drag(r,v);
end
if perturbs.srp
    aXYZ = aXYZ + opts.srp(r,t);
end

% Set the derivatives of the state
dy = [	v
		aGrav+aXYZ ];
end