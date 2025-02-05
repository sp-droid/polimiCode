clear
close all

%% Group requirements
% Nominal orbit, last 3 elements are free to choose
a = 26910; e = 0.7399; i = 62.5711;
% Repeating GT k:m = 2:1
% Perturbation 1: J2
% Perturbation 2: Drag
cD = 2.1;
AoverM = 7.88241304748947 / 900; % m^2/kg

JD = 2.4595195e+06;
jdT = (JD - 2451545.0 ) / 36525;
opts.theta0 = wrapTo2Pi(deg2rad(280.46061837 + 360.98564736629*(JD-2451545.0) ...
    + 0.000387933*jdT^2 - jdT^3/38710000.0));

%% Chosen inputs
bOmega = 74.4297; sOmega = 278.065; theta = 174.3877;
tWindow = [date2mjd2000([2021;11;1;7;19;20.285]);0];
nsteps = 800;
nmax = 360;
cR = 1.25; %Radiation pressure coefficient

%% Inputs
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.09053; % Sidereal day
muSun = astroConstants(4);
Rsun = astroConstants(3);
TTsun = 5778;
muMoon = astroConstants(20);
CS = load('egm96/egm96_to360.ascii', '-ascii');

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
wEarth = 2*pi/Tearth;
[A,B] = legendreAB(nmax);

%% Orbit
% Perform the integration
sunPos = @(t) relativeSun(t, tWindow(1), muSun);
moonPos = @(t) relativeMoon(t, tWindow(1));
opts.wEarth = wEarth;
densityModel = @(r) densitySimplified(norm(r)-Rearth);

opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
%opts.j2Pert = @(r) j2Pert(r,J2,Rearth);
opts.sunThirdBody = @(r,t) thirdBodyPert(r, sunPos(t)-r, muSun);
opts.moonThirdBody = @(r,t) thirdBodyPert(r, moonPos(t)-r, muMoon);
opts.egm96 = @(r,thetaG) egm96(r, thetaG, Rearth, muEarth, nmax, CS, A, B);
opts.relativEffect = @(r,v) relativEffect(r, v, muEarth);
opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, cD, AoverM);
opts.srp = @(r,t) srp(r, sunPos(t)-r, TTsun, Rsun, Rearth, cR, AoverM);
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
    perturbs(j,4) = norm(opts.relativEffect(r,v));
    perturbs(j,5) = norm(opts.drag(r,v));
    perturbs(j,6) = norm(opts.srp(r,t));
end

%% Perturbations
figure;
subplot(3,2,1)
plot(scaledT,perturbs(:,1),'LineWidth',2); title('Sun3rdBody');
subplot(3,2,2)
plot(scaledT,perturbs(:,2),'LineWidth',2); title('Moon3rdBody');
subplot(3,2,3)
plot(scaledT,perturbs(:,3),'LineWidth',2); title('Egm96');
subplot(3,2,4)
plot(scaledT,perturbs(:,4),'LineWidth',2); title('Relativistic');
subplot(3,2,5)
plot(scaledT,perturbs(:,5),'LineWidth',2); title('Drag');
subplot(3,2,6)
plot(scaledT,perturbs(:,6),'LineWidth',2); title('SRP');

%% Static plot
screen = get(0, 'ScreenSize');
[dist,j] = max(vecnorm(Y(:,1:3)'));
r  = normalize([6.8766;-1.7488;0.6759]','norm')*dist; %Get it with ax.CameraPosition
robs = r;
trackT = (rotRz(deg2rad(-angleEarth(j)))*track')';

figure('Color','k','Position', [0 0 screen(3) screen(4)]);
% Celestial bodies
p3Dopts.Units = 'km';
p3Dopts.RotAngle = angleEarth(j);
planet3D('Earth', p3Dopts);
hold on
p3Dopts = rmfield(p3Dopts, 'RotAngle');
p3Dopts.Position = r;
p3Dopts.Size = 3*dist;
planet3D('Milkyway', p3Dopts);
sunRelPos = (normalize(rSun(j,:),'norm')*dist);
sunRelSize = 2.1*Rsun/norm(rSun(j,:)-robs)*norm(sunRelPos-robs)*(1-Rearth/norm(r));
if (dot(robs,sunRelPos)<0)
    [sunX,sunY,sunZ]=sphere;
    sunX = sunX*sunRelSize+sunRelPos(1);
    sunY = sunY*sunRelSize+sunRelPos(2);
    sunZ = sunZ*sunRelSize+sunRelPos(3);
    surf(sunX,sunY,sunZ,'EdgeColor','none','FaceColor',[0.98;0.843;0.627],'AmbientStrength',1,'SpecularStrength',1,'SpecularExponent',500)
end
%Sunlight, apparent Sun is ~28-34 arc minutes near Earth
light("Style","infinite","Position",sunRelPos)
light("Style","local","Position",sunRelPos+normalize(r-sunRelPos,'norm')*2*sunRelSize);
p3Dopts.Position = (normalize(rMoon(j,:),'norm')*dist)';
p3Dopts.Size = norm(p3Dopts.Position-robs)/norm(rMoon(j,:)-robs)*(1-Rearth/norm(r));
if (dot(robs,p3Dopts.Position)<0)
    planet3D('Moon', p3Dopts);
end
% Flight path
scatter3( trackT(:,1), trackT(:,2), trackT(:,3), 12, scaledT, 'filled')
% Sun hitting the s/c
scatter3( Y(:,1), Y(:,2), Y(:,3), 12, scaledT, 'filled')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
%title(datestr(datetime(mjd20002date(tWindow(1)+T(end)/86400))), 'FontSize', 20,'Color','w');
%cbar = colorbar; cbar.Color = 'w'; cbar.Title.Color = 'w';
%cbar.Title.String = strcat('Time [',Tname,']');
clim([min(scaledT);max(scaledT)])   
axis equal;
grid on;
ax = gca; ax.Color = 'k'; ax.GridColor = 'k';
ax.GridAlpha = 0; ax.XColor = 'k'; ax.YColor = 'k';
ax.ZColor = 'k';
ax.CameraPosition = r;  % Set the camera positiont
%camva(2*viewAngle);
ax.XLim = [-dist, dist];  % Set the x-axis range
ax.YLim = [-dist, dist];  % Set the y-axis range
ax.ZLim = [-dist, dist];  % Set the z-axis range
hold off


%% Functions
relativeSun(0,tWindow(1),muSun)
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