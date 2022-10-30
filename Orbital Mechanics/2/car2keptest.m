clc
clear
r = [1.440084852930837e+03;6.464286824897491e+03;-1.518256115919854e+03];
v = [-4.460570434953907;2.344228748475583;5.773283517522664];
mu_E = astroConstants(13);

[a,e,i,bOmega,sOmega,theta] = car2kep(r,v,mu_E,'deg');
% a = 6801.3;
% e = 0.0012;
% i = 0.9038;
% bOmega = 1.5331;
% sOmega = 0.3099;
% theta = 5.6849;
%[r_new,v_new] = kep2car(a,e,i,bOmega,sOmega,theta,mu_E,'deg')

% Function car2kep
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
% VERSIONS
% 2022-10-13: v1
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

% Function kep2car
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

function R1 = rotRx( angle )
%timescaling for relevant time information in plots
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around X axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
R1 = [1 0 0;
    0 cos(angle) sin(angle);
    0 -sin(angle) cos(angle)];
end

function R2 = rotRy( angle )
%timescaling for relevant time information in plots
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around X axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
R2 = [cos(angle) 0 sin(angle);
      0 1 0;
      -sin(angle) 0 cos(angle)];
end

function R3 = rotRz( angle )
%timescaling for relevant time information in plots
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around X axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
R3 = [cos(angle) sin(angle) 0;
      -sin(angle) cos(angle) 0;
      0 0 1];
end