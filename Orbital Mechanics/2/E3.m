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
img = imread('earth2D','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on

[long,lat] = removeLonLines(long,lat);
plot(long,lat,'green')
hold on

%% Propagating from initial ephemerids
Torb = 2*pi*sqrt( isseph(1,7)^3/mu_E );

nOrb = 60*length(isseph)/Torb;
nPoints = 4000;

y0 = [r(1,:) v(1,:)];

% Perform the integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[ Y, t ] = timed2BP(y0,mu_E,opts,nPoints,60*length(isseph));
t = t';

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