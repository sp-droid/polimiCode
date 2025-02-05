% Evolution of keplerian elements
clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

%% Inputs

% Celestial parameters
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.1;
wEarth = 2*pi/Tearth;
J2 = astroConstants(9);
greenwich = 0;

muSun = astroConstants(4);
Rsun = astroConstants(3);
TTsun = 5778;
muMoon = astroConstants(20);

cD = 2.7;
AoverMdrag = 0.0045; % m^2/kg
AoverMsrp = 0.016;

% Orbit parameters
a = 7271;                           % Semi-major axis
e = 0.001;                          % Eccentricity
i = 78;                         % Inclination
bOmega = 0;                         % Right ascension of ascending node
sOmega = 0;                         % Argument of pericentre
theta = 0;                           % Initial true anomaly

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];

date = [2032;1;1;0;0;0];
tWindow = [date2mjd2000(date);0];
norbits = 10;
ngrid = 100000;

%% Orbit propagation
sunPos = @(t) relativeSun(t, tWindow(1), muSun);
moonPos = @(t) relativeMoon(t, tWindow(1));
opts.wEarth = wEarth;
cR = 1.4; % Reflectivity coefficient
Ap = 10;  % Geomagnetic index

CS = load('egm96_to360.ascii', '-ascii');
nmax = 360;
[A,B] = legendreAB(nmax);
JD = date2jd(date);
jdT = (JD - 2451545.0 ) / 36525;
opts.theta0 = wrapTo2Pi(deg2rad(280.46061837 + 360.98564736629*(JD-2451545.0) ...
    + 0.000387933*jdT^2 - jdT^3/38710000.0));

opts.sunThirdBody = @(r,t) thirdBodyPert(r, sunPos(t)-r, muSun);
%opts.moonThirdBody = @(r,t) thirdBodyPert(r, moonPos(t)-r, muMoon);
%opts.j2Pert = @(r) j2Pert(r,0.0010826,Rearth,muEarth);
%opts.egm96 = @(r,thetaG) egm96(r, thetaG, Rearth, muEarth, nmax, CS, A, B);
%opts.relativEffect = @(r,v) relativEffect(r, v, muEarth);
%opts.drag = @(r,v,t) drag(r, v, densityHigh(norm(r)-Rearth,f107estimation(t,tWindow(1)),Ap), wEarth, cD, AoverMdrag);
%opts.srp = @(r,t) srp(r, sunPos(t)-r, TTsun, Rsun, Rearth, cR, AoverMsrp);
%opts.TinPeriods = true;
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[Y1car, T] = timed2BP(y0, muEarth, opts, ngrid, [0, 63072000]);

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

[scaledT, Tname] = timescaling(T);
%% Plots

figure;
subplot(2,3,1)
triPlot(T, Y1kep(:,1), 'a axis [km]', norbits)
subplot(2,3,2)
triPlot(T, Y1kep(:,2), 'e [-]', norbits)
subplot(2,3,3)
triPlot(T, rad2deg(Y1kep(:,3)), 'i [deg]', norbits)
subplot(2,3,4)
triPlot(T, rad2deg(Y1kep(:,4)), 'Ω [deg]', norbits, true)
subplot(2,3,5)
triPlot(T, rad2deg(Y1kep(:,5)), 'ω [deg]', norbits, true)
subplot(2,3,6)
triPlot(T, rad2deg(Y1kep(:,6)), 'ϑ [deg]', norbits, true)
print(gcf, 'output/pertSun.png', '-dpng', '-r300');

% Plot percentage difference in 100T, value in 100T and value in first 10%
function triPlot(T, var, name, norbs, xlabelBool)
if nargin < 5
    xlabelBool = false; % Default value for opt1
end

[scaledT, Tname] = timescaling(T);

% nsecular = ceil(1*length(var)/norbs); % 1 orbits
% nlongterm = ceil(0.25*length(var)/norbs); % 0.25 orbits
% nshortterm = ceil(0.05*length(var)/norbs); % 0.05 orbits
% 
% secular = movmean( var, nsecular );
movWindow = 100;
longterm = movmean( var, movWindow );
% shortterm = movmean( var, nshortterm );

plot(scaledT, var,'LineWidth',1)
hold on
plot(scaledT, longterm,'LineWidth',1)
%plot(T(nlongterm:end-nlongterm), longterm(nlongterm:end-nlongterm),'LineWidth',1)
hold on
%plot(scaledT(nsecular:end-nsecular), secular(nsecular:end-nsecular),'LineWidth',1)
% plot(T(nshortterm:end-nshortterm), shortterm(nshortterm:end-nshortterm),'LineWidth',0.5)

if (xlabelBool == true)
	xlabel(strcat('time [',Tname,']'));
end
title(name);
if (strcmp(name,'True anomaly [deg]'))
    legend('Value',strcat('Moving window n=',string(movWindow)))
end
grid on;
axis tight;
hold off

% secular = secular(nsecular:end-nsecular);
% X = [ones(length(secular),1) T(nsecular:end-nsecular)];
% b = X\secular;
% R2 = 1-sum((secular - X*b).^2)/sum((secular - mean(secular)).^2);
%fprintf(strcat(name,": ",string(b(1))," ",string(b(2)),"*orbit, R2: ",string(R2),"\n"));
end

function r = relativeSun(tseconds, initialMJD, muSun)
    T = initialMJD + tseconds/24/3600;
    [kep,~] = uplanet(T, 3);
    [r,~] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
    % This is a position vector pointing from the Sun to the Earth, but in
    % ecliptic frame. We need to transform back to equatorial
    r = -r;
    
    tilt = deg2rad(23.4365472133); %Ecliptic obliquity in January 2021
    r = rotRx(tilt)'*r;
end

function r = relativeMoon(tseconds, initialMJD)
    T = initialMJD + tseconds/24/3600;
    [r,~] = ephMoon(T);
    r = r'; % 3x1
end