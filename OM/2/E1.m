clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);
% Satellite groundtrack

% Physical parameters
mu_E = astroConstants(13);
R_E = astroConstants(23);
w_E = deg2rad(15.04/3600);
greenwich = 0;

% Orbit parameters
r0 = [ -4578.219; -801.084; -7929.708]; % [km]
v0 = [ 0.8; -6.037; 1.385]; % [km/s]
y0 = [ r0; v0 ];

a = 8350;                           % Semi-major axis
e = 0.1976;                         % Eccentricity
i = deg2rad(60);                    % Inclination
bOmega = deg2rad(270);              % Right ascension of ascending node
sOmega = deg2rad(45);               % Argument of pericentre
f0 = deg2rad(230);                  % Initial true anomaly

Torb = 2*pi*sqrt( a^3/mu_E );          % Orbital period
lambda = Torb*w_E;                     % Ground track drift
nOrb = 3.25;
nPoints = 5000;
t = linspace(0,Torb*nOrb,nPoints);

% NOTA: IMPLEMENTAR FUNCION QUE TRANSFORME STATE -> KEP VECTOR
%% Computation
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E, 0, R_E), t, y0, options );

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');
nOrb = t/Torb;


delta = asin(r(:,3)'./Rnorm);                  % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-w_E*t));
lat = rad2deg(delta);

%% Plots
figure()
img = imread('earth2D','png');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on
scatter(long, lat, 2, nOrb)
xlim([-180,180]);
ylim([-90,90]);
colorbar;
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