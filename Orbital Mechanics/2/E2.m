clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);
% Satellite with REPEATING groundtrack

% Repeating conditions
k = 12;                             % Number of satellite orbits
m = 1;                              % Number of planetary rotations

% Physical parameters
mu_E = astroConstants(13);
R_E = astroConstants(23);
w_E = deg2rad(15.04/3600);
greenwich = 0;

% Orbit parameters
r0 = [ -4578.219; -801.084; -7929.708]; % [km]
v0 = [ 0.8; -6.037; 1.385]; % [km/s]
y0 = [ r0; v0 ];

a = 8350;                           % Semi-major axis of UNMODIFIED ORBIT
e = 0.1976;                         % Eccentricity
i = deg2rad(60);                    % Inclination
bOmega = deg2rad(270);              % Right ascension of ascending node
sOmega = deg2rad(45);               % Argument of pericentre
f0 = deg2rad(230);                  % Initial true anomaly

Torb = 2*pi*sqrt( a^3/mu_E );       % Orbital period
lambda = Torb*w_E;                  % Ground track drift

% Repeating ground track condition
n = w_E*k/m;                        % Mean motion of modified orbit
a_RGT = (mu_E/n^2)^(1/3);           % Modified semi-major axis

Torb_RGT = 2*pi*sqrt( a_RGT^3/mu_E );

nOrb = k;
nPoints = 15000;

%% Computation
% Perform the integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[ Y, t ] = timed2BP(y0,mu_E,opts,nPoints,Torb_RGT*nOrb);
t = t';

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');
nOrb = t/Torb;


delta = asin(r(:,3)'./Rnorm);                   % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-w_E*t));
lat = rad2deg(delta);

figure()
img = imread('earth2D','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on

scatter(long, lat, 3, 'green')
hold on

%% Repeating groundtrack
% Let's say we go from one orbit to the other by modifying V at the same
% starting position, so r remains the same.
fact = sqrt(mu_E*(2/vecnorm(r0)-1/a_RGT))/vecnorm(v0);
v0_RGT = v0.*fact;
y0_RGT = [ r0; v0_RGT ];

% Perform the integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[ Y, t ] = timed2BP(y0_RGT,mu_E,opts,nPoints,[],k);
t = t';


r_RGT = Y(:,1:3);
v_RGT = Y(:,4:6);
Rnorm_RGT = vecnorm(r_RGT');
Vnorm_RGT = vecnorm(v_RGT');
nOrb_RGT = t/Torb_RGT;

delta_RGT = asin(r_RGT(:,3)'./Rnorm_RGT);                   % Declination
alpha_RGT = atan2(r_RGT(:,2)',r_RGT(:,1)');                 % Right ascension
long_RGT = wrapTo180(rad2deg(alpha_RGT-greenwich-w_E*t));
lat_RGT = rad2deg(delta_RGT);

scatter(long_RGT, lat_RGT, 3, 'red')
hold on

plot(long(1),lat(1),'^','Color',[0,1,0],'LineWidth',6)
plot(long(end),lat(end),'v','Color',[0,0.8,0],'LineWidth',6)
plot(long_RGT(1),lat_RGT(1),'^','Color',[1,0,0],'LineWidth',4)
plot(long_RGT(end),lat_RGT(end),'v','Color',[0.8,0,0],'LineWidth',4)

legend('Original [OG]','Repeating [RGT]','OG start','OG finish','RGT start', 'RGT finish','Location','northoutside','NumColumns',6)
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
grid on
hold off