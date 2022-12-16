function [r,v] = kep2car( a, e, i, bOmega, sOmega, theta, mu, angleUnit )
%Conversion from cartesian coordinates to keplerian elements
%
% PROTOTYPE
% r, v = car2kep( a, e, i, bOmega, sOmega, theta, mu_E, angleUnit )
%
% INPUT:
% a[1] Semi-major axis
% e[1] Eccentricity
% i[1] Inclination
% bOmega[1] Right ascension of the ascending node
% sOmega[1] Argument of periapsis
% theta[1] True anomaly
% mu[1] Standard gravitational parameter
% angleUnit[str] Possibles 'rad' or 'deg'. Radians by default
%
% OUTPUT:
% r[3x1] Position vector
% v[3x1] Velocity vector
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
%Conversion to other angle units if it's necessary
if isequal(angleUnit,'deg')
    i = deg2rad(i);
    bOmega = deg2rad(bOmega);
    sOmega = deg2rad(sOmega);
    theta = deg2rad(theta);
end
p = a*(1-e^2);
hNorm = sqrt(mu*p);
rNorm = p/(1+e*cos(theta));

%r and v in perifocal frame {e, p, h}
%r = rNorm*[cos(theta); sin(theta); 0];
%v = mu/hNorm*[-sin(theta); e+cos(theta); 0];
%r and v in orbital frame {r, theta, h}
r = [rNorm; 0; 0];
v = mu/hNorm*[e*sin(theta); 1+e*cos(theta); 0];

%Transformation of perifocal coordinates into inertial equatorial ones
%A = (rotz(sOmega)*rotx(i)*rotz(bOmega))';
%Transformation of orbital coordinates into inertial equatorial ones
A = (rotRz(sOmega+theta)*rotRx(i)*rotRz(bOmega))';

r = A*r;
v = A*v;
end