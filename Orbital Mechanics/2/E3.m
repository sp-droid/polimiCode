clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);
% ISS groundtrack

% Physical parameters
mu_E = astroConstants(13);
R_E = astroConstants(23);
w_E = deg2rad(15.04/3600);
greenwich = 0;

% The raw iss_eph text file was filtered to a 7 column matrix already
% Load ephemerids (already saved as matlab file with "save filename var")
load("iss_ephemerids.mat")

%% Plotting the downloaded ephemerids
t = linspace(0,60*length(isseph),length(isseph));
r = zeros(length(isseph),3);
v = zeros(length(isseph),3);
for j=1:length(isseph)
    [r(j,:),v(j,:)] = kep2car(isseph(j,7),isseph(j,2),isseph(j,3),isseph(j,4),isseph(j,5),isseph(j,6),mu_E,'deg');
end

Rnorm = vecnorm(r');
Vnorm = vecnorm(v');

delta = asin(r(:,3)'./Rnorm);                   % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-w_E*t))';
lat = rad2deg(delta)';

figure()
img = imread('earth2D','png');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on

[long,lat] = removeLonLines(long,lat);
plot(long,lat,'green')
hold on

%% Propagating from initial ephemerids
Torb = 2*pi*sqrt( isseph(1,7)^3/mu_E );

nOrb = 60*length(isseph)/Torb;
nPoints = 4000;
t = linspace(0,60*length(isseph),nPoints);

y0 = [r(1,:) v(1,:)];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E, 0, R_E), t, y0, options );

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');

delta = asin(r(:,3)'./Rnorm);                   % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long_prop = wrapTo180(rad2deg(alpha-greenwich-w_E*t))';
lat_prop = rad2deg(delta)';

[long_prop,lat_prop] = removeLonLines(long_prop,lat_prop);
plot(long_prop, lat_prop, 'red','LineWidth',2)
hold on

plot(long(1),lat(1),'^','Color',[0,1,0],'LineWidth',6)
plot(long(end),lat(end),'v','Color',[0,0.8,0],'LineWidth',6)
plot(long_prop(1),lat_prop(1),'^','Color',[1,0,0],'LineWidth',4)
plot(long_prop(end),lat_prop(end),'v','Color',[0.8,0,0],'LineWidth',4)

legend('Ephemerids','Propagated','Eph. start','Eph. finish','Prop. start', 'Prop. finish','Location','northoutside','NumColumns',6)
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on
hold off

%% Functions
% Function ode_2bp
function dy = ode_2bp( ~, y, mu, J2, R )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);

% aJ2 term
aJ2 = 5*r(3)^2/rnorm^2;
aJ2 = [r(1)/rnorm*(aJ2-1)
    r(2)/rnorm*(aJ2-1)
    r(3)/rnorm*(aJ2-3)];
aJ2 = aJ2*1.5*J2*mu*R^2/rnorm^4;

% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r+aJ2 ];
end

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
%Returns a rotation matrix around the first axis
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
%Returns a rotation matrix around the second axis
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
%Returns a rotation matrix around the third axis
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

function [newlong,newlat] = removeLonLines( long, lat )
%Adds NaN values every time longitude wraps around the planet for better line plots
%
% PROTOTYPE
% long = rotx( long )
%
% INPUT:
% long[nx1] Longitude in degrees
% lat[nx1] Latitude in any unit
%
% OUTPUT:
% newlong[(n+k)x1] Longitude with intermediate NaN values
% newlat[nx1] Latitude with intermediate NaN values
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-14: v1
%
% -------------------------------------------------------------------------
newlong = [];
newlat = [];
for j=1:length(long)-1
    if long(j)>150 && long(j+1)<30
        newlong = [newlong; long(j); NaN];
        newlat = [newlat; lat(j); NaN];
    elseif long(j)<30 && long(j+1)>150
        newlong = [newlong; long(j); NaN];
        newlat = [newlat; lat(j); NaN];
    else
        newlong = [newlong; long(j)];
        newlat = [newlat; lat(j)];
    end
end
newlong = [newlong; long(end)];
newlat = [newlat; lat(j)];
end