clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Inputs
% Celestial parameters
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

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
tWindow = [date2mjd2000([2021;11;1;0;0;0]);0];
norbits = 10;
ngrid = 10000;

[a,e,i,bOmega,sOmega,theta] = car2kep(r0,v0,muEarth,'rad');
kep0 = [a,e,i,bOmega,sOmega,theta]';

%% Orbit propagation
densityModel = @(r) densitySimplified(norm(r)-Rearth);
opts.j2Pert = @(r) j2Pert(r,J2,Rearth,muEarth);
opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, 2.1, 0.0095);
opts.TinPeriods = true;
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[Y1car, T] = timed2BP(y0, muEarth, opts, ngrid, [], norbits);
opts.keplerian = true;
[Y2kep, ~] = timed2BP(kep0, muEarth, opts, ngrid, [], norbits);
Y1kep = zeros(ngrid,6);
for j=1:ngrid
    [a,e,i,bOmega,sOmega,theta] = car2kep(Y1car(j,1:3),Y1car(j,4:6),muEarth,'rad');
    Y1kep(j,1) = a; Y1kep(j,2) = e; Y1kep(j,3) = i;
    Y1kep(j,4) = bOmega; Y1kep(j,5) = sOmega; Y1kep(j,6) = theta;
end

% Unwrap keplerian 4, 6 and 6
Y1kep(:,4) = unwrap(Y1kep(:,4));
Y1kep(:,5) = unwrap(Y1kep(:,5));
Y1kep(:,6) = unwrap(Y1kep(:,6));
Y2kep(:,4) = unwrap(Y2kep(:,4));
Y2kep(:,5) = unwrap(Y2kep(:,5));
Y2kep(:,6) = unwrap(Y2kep(:,6));

[scaledT, Tname] = timescaling(T);
%% Plots
figure;
subplot(3,2,1)
plot(scaledT,(Y1kep(:,1)-Y2kep(:,1))/a,'LineWidth',2); title('Semi-major axis');
subplot(3,2,2)
plot(scaledT,(Y1kep(:,2)-Y2kep(:,2))/1,'LineWidth',2); title('Eccentricity');
subplot(3,2,3)
plot(scaledT,(Y1kep(:,3)-Y2kep(:,3))/pi,'LineWidth',2); title('Inclination');
subplot(3,2,4)
plot(scaledT,(Y1kep(:,4)-Y2kep(:,4))/(2*pi),'LineWidth',2); title('RAAN');
subplot(3,2,5)
plot(scaledT,(Y1kep(:,5)-Y2kep(:,5))/(2*pi),'LineWidth',2); title('Argument of periapsis');
subplot(3,2,6)
plot(scaledT,(Y1kep(:,6)-Y2kep(:,6))/(2*pi),'LineWidth',2); title('True anomaly');

figure;
subplot(3,2,1)
triPlot(T, Y1kep(:,1), 'Semi-major axis', norbits)
subplot(3,2,2)
triPlot(T, Y1kep(:,2), 'Eccentricity', norbits)
subplot(3,2,3)
triPlot(T, rad2deg(Y1kep(:,3)), 'Inclination', norbits)
subplot(3,2,4)
triPlot(T, rad2deg(Y1kep(:,4)), 'RAAN', norbits)
subplot(3,2,5)
triPlot(T, rad2deg(Y1kep(:,5)), 'Argument of periapsis', norbits)
subplot(3,2,6)
triPlot(T, rad2deg(Y1kep(:,6)), 'True anomaly', norbits)

% Plot percentage difference in 100T, value in 100T and value in first 10%
function triPlot(T, var, name, norbs)

secular = movmean( var, ceil(length(var)/norbs) ); % 1 orbit
longterm = movmean( var, ceil(0.25*length(var)/norbs) ); % 0.25 orbits
shortterm = movmean( var, ceil(0.05*length(var)/norbs) ); % 0.05 orbits

%plot(T, var)
hold on
plot(T, secular,'LineWidth',2)
hold on
plot(T, longterm,'LineWidth',2)
plot(T, shortterm)
xlabel('time [T]'); ylabel('var value');
title(name);
legend('Secular','Long term','Short term')
grid on;
axis tight;
hold off

name
X = [ones(length(T),1) T];
b = X\secular
R2 = 1-sum((secular - X*b).^2)/sum((secular - mean(secular)).^2)
end