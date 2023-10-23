clear
close all
%set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);
% Satellite groundtrack

% Physical parameters
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.1;
wEarth = 2*pi/Tearth;
J2 = astroConstants(9);
greenwich = 0;

% Orbit parameters
a = 26619;                           % Semi-major axis
e = 0.7452;                          % Eccentricity
i = 62.9089;                         % Inclination
bOmega = 60;                         % Right ascension of ascending node
sOmega = 30;                         % Argument of pericentre
theta = 0;                           % Initial true anomaly

% Repeating GT
koverm = 2;
a = ((Tearth/koverm/2/pi)^2*muEarth)^(1/3);

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
Torb = 2*pi*sqrt( a^3/muEarth );          % Orbital period
lambda = Torb*wEarth;                     % Ground track drift
nOrb = 1;
nPoints = 35000;
ttime = 20*Tearth;

%% Computation
% Perform the integration
opts.RelTol = 1e-12;
opts.AbsTol = 1e-13;

[ Y, t ] = timed2BP(y0,muEarth,opts,nPoints,ttime);
t = t';

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');
tOrbs = t/Torb;


delta = asin(r(:,3)'./Rnorm);                  % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-wEarth*t));
lat = rad2deg(delta);

%% Plots
figure('Position', [0 0 1080 270]);
subplot(1,2,1)
img = imread('earth2D','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on
scatter(long, lat, 3, 'MarkerEdgeColor','red','MarkerFaceColor','red')%tOrbs)
hold on

plot(long(1),lat(1),'^','Color',[0,1,0],'LineWidth',4)
plot(long(end),lat(end),'v','Color',[0,0,1],'LineWidth',4)

%legend('Groundtrack','Start','Finish','Location','northoutside','NumColumns',3)
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
%hcb=colorbar;
%title(hcb,'Orbit number')
grid on
hold off

%%
densityModel = @(r) densitySimplified(norm(r)-Rearth);
opts.j2Pert = @(r) j2Pert(r,J2,Rearth,muEarth);
opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, 2.1, 0.0095);

[ Y, t ] = timed2BP(y0,muEarth,opts,nPoints,ttime);
t = t';

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');
tOrbs = t/Torb;


delta = asin(r(:,3)'./Rnorm);                  % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-wEarth*t));
lat = rad2deg(delta);

subplot(1,2,2)
img = imread('earth2D','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on
scatter(long, lat, 3, 'MarkerEdgeColor','red','MarkerFaceColor','red')%tOrbs)
hold on

plot(long(1),lat(1),'^','Color',[0,1,0],'LineWidth',4)
plot(long(end),lat(end),'v','Color',[0,0,1],'LineWidth',4)

%legend('Groundtrack','Start','Finish','Location','northoutside','NumColumns',3)
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
%hcb=colorbar;
%title(hcb,'Orbit number')
grid on
hold off