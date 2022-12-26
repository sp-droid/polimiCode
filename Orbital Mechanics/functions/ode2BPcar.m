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
if perturbs.J2
	aJ2 = 5*r(3)^2/rnorm^2;
	aJ2 = [	r(1)/rnorm*(aJ2-1)
			r(2)/rnorm*(aJ2-1)
			r(3)/rnorm*(aJ2-3)];
	aJ2 = aJ2*1.5*opts.J2*mu*opts.Rearth^2/rnorm^4;
	aXYZ = aXYZ + aJ2;
end
if perturbs.sun
    rSun = opts.sunPos(t);
    aSun = opts.muSun*((rSun-r)/norm(rSun-r)^3-rSun/norm(rSun)^3);
    aXYZ = aXYZ + aSun;
end
if perturbs.moon
    rMoon = opts.moonPos(t);
    aMoon = opts.muMoon*((rMoon-r)/norm(rMoon-r)^3-rMoon/norm(rMoon)^3);
    aXYZ = aXYZ + aMoon;
end
if perturbs.egm96
    % Correction
    thetaG = opts.wEarth*(t-opts.t0);
    aXYZ = aXYZ + opts.egm96(r, thetaG);
end
if perturbs.relativ
    c = opts.lightSpeed;
    aGrav = aGrav - mu/rnorm^3*(((2/c)^2*mu/rnorm-dot(v,v)/c^2)*r+(2/c)^2*dot(r,v)*v);
end

% Set the derivatives of the state
dy = [	v
		aGrav+aXYZ ];
end