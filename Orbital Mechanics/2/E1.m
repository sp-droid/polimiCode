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
%r0 = [ -4578.219; -801.084; -7929.708]; % [km]
%v0 = [ 0.8; -6.037; 1.385]; % [km/s]
r0 = [ 1599.4; 5859.1; 3019.2 ]; % [km]
v0 = [ -5.9909; -2.3882; 7.8083 ]; % [km/s]
y0 = [ r0; v0 ];

% a = 8350;                           % Semi-major axis
% e = 0.1976;                         % Eccentricity
% i = deg2rad(60);                    % Inclination
% bOmega = deg2rad(270);              % Right ascension of ascending node
% sOmega = deg2rad(45);               % Argument of pericentre
% f0 = deg2rad(230);                  % Initial true anomaly
[a,e,i,bOmega,sOmega,f0] = car2kep(r0,v0,mu_E,'rad');

Torb = 2*pi*sqrt( a^3/mu_E );          % Orbital period
lambda = Torb*w_E;                     % Ground track drift
nOrb = 5.5;
nPoints = 5000;

% NOTA: IMPLEMENTAR FUNCION QUE TRANSFORME STATE -> KEP VECTOR
%% Computation
% Perform the integration
opts.RelTol = 1e-12;
opts.AbsTol = 1e-13;
[ Y, t ] = timed2BP(y0,mu_E,opts,nPoints,Torb*nOrb);
t = t';

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
img = imread('earth2D','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on
scatter(long, lat, 2, nOrb)
hold on

plot(long(1),lat(1),'^','Color',[0,1,0],'LineWidth',6)
plot(long(end),lat(end),'v','Color',[0,0.8,0],'LineWidth',6)

legend('Groundtrack','Start','Finish','Location','northoutside','NumColumns',3)
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
hcb=colorbar;
title(hcb,'Orbit number')
grid on
hold off