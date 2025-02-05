clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

%% Group requirements
% Nominal orbit, last 3 elements are free to choose
a = 7271; e = 0.0; i = 78;
% Repeating GT k:m = 2:1
% Perturbation 1: J2
% Perturbation 2: Drag
cD = 2.7;
AoverMdrag = 0.0045; % m^2/kg
AoverMsrp = 0.016;

%% Chosen inputs
bOmega = 0; sOmega = 0; theta = 0;
date = [2032;1;1;0;0;0];
tWindow = [date2mjd2000(date);0];
nsteps = 200;
nmax = 360;
cR = 1.4; %Radiation pressure coefficient
Ap = 10;

%% Inputs
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.09053; % Sidereal day
muSun = astroConstants(4);
Rsun = astroConstants(3);
TTsun = 5778;
muMoon = astroConstants(20);
CS = load('egm96_to360.ascii', '-ascii');

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
wEarth = 2*pi/Tearth;
[A,B] = legendreAB(nmax);
JD = date2jd(date);
jdT = (JD - 2451545.0 ) / 36525;
opts.theta0 = wrapTo2Pi(deg2rad(280.46061837 + 360.98564736629*(JD-2451545.0) ...
    + 0.000387933*jdT^2 - jdT^3/38710000.0));

%% Orbit
% Perform the integration
sunPos = @(t) relativeSun(t, tWindow(1), muSun);
moonPos = @(t) relativeMoon(t, tWindow(1));
opts.wEarth = wEarth;

opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
opts.sunThirdBody = @(r,t) thirdBodyPert(r, sunPos(t)-r, muSun);
opts.moonThirdBody = @(r,t) thirdBodyPert(r, moonPos(t)-r, muMoon);
opts.j2Pert = @(r) j2Pert(r,1082.64e-6,Rearth,muEarth);
opts.egm96 = @(r,thetaG) egm96(r, thetaG, Rearth, muEarth, nmax, CS, A, B);
opts.relativEffect = @(r,v) relativEffect(r, v, muEarth);
opts.drag = @(r,v,t) drag(r, v, densityHigh(norm(r)-Rearth,f107estimation(t,tWindow(1)),Ap), wEarth, cD, AoverMdrag);
opts.srp = @(r,t) srp(r, sunPos(t)-r, TTsun, Rsun, Rearth, cR, AoverMsrp);
opts.perturbShow = true;
[ Y, T ] = timed2BP(y0,muEarth,opts,nsteps,[],[0,2]);

% Define time scale, time window and time in mjd
[scaledT, Tname] = timescaling(T);

tWindow(2) = tWindow(1)+T(end)/24/3600;
Tmjd = linspace(tWindow(1), tWindow(2), nsteps)';

%% Celestial bodies
rSun = zeros(nsteps,3);
rMoon = zeros(nsteps,3);
angleEarth = zeros(nsteps,1);
track = zeros(nsteps,3); % This is the track relative to the initial Earth.
% To obtain the one corresponding to time t(j), rotate the Earth angleEarth(j) and
% trackT = (rotRz(deg2rad(-angleEarth(j)))*track')';
% Obtain perturbation numbers for every step
perturbs = zeros(nsteps,6);
for j=1:nsteps
    r = Y(j,1:3)'; v = Y(j,4:6)'; t = T(j); thetaG = wEarth*(t-0);
    rSun(j,:) = relativeSun(t, tWindow(1), muSun);
    angleEarth(j) = rad2deg(wrapTo2Pi(thetaG));
    track(j,1:3) = (rotRz(deg2rad(angleEarth(j)))*normalize(r,'norm'))'*Rearth;
    rMoon(j,:) = relativeMoon(t, tWindow(1));
    perturbs(j,1) = norm(opts.sunThirdBody(r,t));
    perturbs(j,2) = norm(opts.moonThirdBody(r,t));
    perturbs(j,3) = norm(opts.egm96(r,thetaG));
    perturbs(j,4) = norm(opts.egm96(r,thetaG)-opts.j2Pert(r));
    perturbs(j,5) = norm(opts.drag(r,v,t));
    perturbs(j,6) = norm(opts.srp(r,t));
end
perturbs = perturbs*1000;

%% Perturbations
figure;
subplot(2,3,1)
plot(scaledT,perturbs(:,1),'LineWidth',2); title('Sun3rdBody'); grid on; xlim([0,max(scaledT)]);
subplot(2,3,2)
plot(scaledT,perturbs(:,2),'LineWidth',2); title('Moon3rdBody'); grid on; xlim([0,max(scaledT)]);
subplot(2,3,3)
plot(scaledT,perturbs(:,3),'LineWidth',2); title('EGM96'); grid on; xlim([0,max(scaledT)]);
subplot(2,3,4)
plot(scaledT,perturbs(:,4),'LineWidth',2); title('EGM96-J2'); grid on; xlabel(Tname); xlim([0,max(scaledT)]);
subplot(2,3,5)
plot(scaledT,perturbs(:,5),'LineWidth',2); title('Drag'); grid on; xlabel(Tname); xlim([0,max(scaledT)]);
subplot(2,3,6)
plot(scaledT,perturbs(:,6),'LineWidth',2); title('SRP'); grid on; xlabel(Tname); xlim([0,max(scaledT)]);

print(gcf, 'output/perturbations.png', '-dpng', '-r300');

%% Functions
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