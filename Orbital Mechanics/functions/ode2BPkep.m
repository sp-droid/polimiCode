function dy = ode2BPkep( ~, y, mu, opts, perturbs )
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

% Primary gravity acceleration
aGrav = h/rnorm^2;

aRSW = [0;0;0];
% J2 perturbation directly in RSW frame
if perturbs.J2
	aJ2 = [ 1-3*sin(i)^2*sin(theta+sOmega)^2
			sin(i)^2*sin(2*(theta+sOmega))
			sin(2*i)*sin(theta+sOmega)];
	aJ2 = aJ2 * -1.5*J2*mu*Rearth^2/rnorm^4;
	aRSW = aRSW + aJ2;
end

% Set the derivatives of the state in RSW frame
dy = [	2*a^2/h * (e*sin(theta)*aRSW(1) + p/rnorm*aRSW(2))
		1/h * (p*sin(theta)*aRSW(1) + ((p+rnorm)*cos(theta) + rnorm*e)*aRSW(2))
		rnorm/h*cos(theta+sOmega)*aRSW(3)
		rnorm/h*sin(theta+sOmega)/sin(i)*aRSW(3)
		1/h/e * (-p*cos(theta)*aRSW(1) + (p+rnorm)*sin(theta)*aRSW(2)) - rnorm/h*sin(theta+sOmega)*cos(i)/sin(i)*aRSW(3)
		aGrav + 1/h/e * (p*cos(theta)*aRSW(1) - (p+rnorm)*sin(theta)*aRSW(2))];
end