function [a,e,i,bOmega,sOmega,theta] = car2kep( r, v, mu, angleUnit )
%Conversion from cartesian coordinates to keplerian elements
%
% PROTOTYPE
% a, e, i, bOmega, sOmega, theta = car2kep( r, v, mu_E, 'deg' )
%
% INPUT:
% r[3x1] Position vector
% v[3x1] Velocity vector
% mu[1] Standard gravitational parameter
% angleUnit[str] Possibles 'rad' or 'deg'. Radians by default
%
% OUTPUT:
% a[1] Semi-major axis
% e[1] Eccentricity
% i[1] Inclination
% bOmega[1] Right ascension of the ascending node
% sOmega[1] Argument of periapsis
% theta[1] True anomaly
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------
rNorm = vecnorm(r);
vNorm = vecnorm(v);
rVersor = r/rNorm;
a = mu/(2*mu/rNorm-vNorm^2);

h = cross(r,v);

i = acos(h(3)/vecnorm(h));

eVersor = 1/mu*cross(v,h)-rVersor;
e = vecnorm(eVersor);

N = cross([0;0;1],h);
N = N/vecnorm(N);
%If inclination is 0, N / line of nodes is not defined.
% Convention -> N = [1;0;0]
if abs(i) < 1e-7
    N = [1;0;0];
end

%In circular orbits, e / line of apses is not defined.
% Convention -> e = N
if e < 1e-9
    eVersor = N;
else
    eVersor = eVersor/e;
end

theta = acos(dot(eVersor,rVersor));
if dot(r,v) < 0
    theta = 2*pi-theta;
end

bOmega = acos(N(1));
if N(2)<0
    bOmega = 2*pi-bOmega;
end
sOmega = acos(dot(N,eVersor));
if eVersor(3)<0
    sOmega = 2*pi-sOmega;
end

%Conversion to other angle units if it's necessary
if isequal(angleUnit,'deg')
    i = rad2deg(i);
    bOmega = rad2deg(bOmega);
    sOmega = rad2deg(sOmega);
    theta = rad2deg(theta);
end
end