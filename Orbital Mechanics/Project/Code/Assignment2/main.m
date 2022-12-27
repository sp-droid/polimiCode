clear
close all

%% Group requirements
% Nominal orbit, last 3 elements are free to choose
a = 26619; e = 0.7452; i = 62.9089;
% Repeating GT k:m = 2:1
% Perturbation 1: J2
% Perturbation 2: Drag
cD = 2.1;
AoverM = 0.0095; % m^2/kg

%% Chosen inputs
bOmega = 60; sOmega = 30; theta = 0;
tWindow = [date2mjd2000([2023;11;1;0;0;0]);0];
nsteps = 2400;
nmax = 360;
cR = 1; %Radiation pressure coefficient

%% Inputs
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.1; % Sidereal day
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
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
%opts.sunPos = @(t) relativeSun(t, tWindow(1), muSun);
%opts.muSun = muSun;
%opts.moonPos = @(t) relativeMoon(t, tWindow(1));
%opts.muMoon = muMoon;
%opts.wEarth = wEarth;
%opts.egm96 = @(r, thetaG) egm96acc(r, thetaG, Rearth, muEarth, nmax, CS, A, B);
%opts.relativ = true;
%densityModel = @(r) densitySimplified(norm(r)-Rearth);
%opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, cD, AoverM);
%sunPos = @(t) relativeSun(t, tWindow(1), muSun);
%opts.srp = @(r,t) srp(r, sunPos, TTsun, Rsun, cR, AoverM)
%i should clean this up a bit, make sure every perturbation is 1 function inside timed2bp
opts.perturbShow = true;
[ Y, T ] = timed2BP(y0,muEarth,opts,nsteps,[],1);

tt = vecnorm(Y(:,1:3)'); tt(end);
tt = max(vecnorm(Y(:,4:6)')); tt(end);
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
sunP = zeros(nsteps,1);
for j=1:nsteps
    rSun(j,:) = relativeSun(T(j), tWindow(1), muSun);
    angleEarth(j) = rad2deg(wrapTo2Pi(wEarth*T(j)));
    track(j,1:3) = (rotRz(deg2rad(angleEarth(j)))*normalize(Y(j,1:3)','norm'))'*Rearth;
    rMoon(j,:) = relativeMoon(T(j), tWindow(1));
    sunP(j) = norm(srp(Y(j,1:3),rSun(j,:),TTsun,Rsun,Rearth,cR,AoverM));
end

%% Static plot
screen = get(0, 'ScreenSize');
j = nsteps;
r  = Y(j,1:3);
trackT = (rotRz(deg2rad(-angleEarth(j)))*track')';
dist = norm(r);%tand(68)*(norm(r)-Rearth);
%
figure('Color','k','Position', [0 0 screen(3) screen(4)]);
% Celestial bodies
p3Dopts.Units = 'km';
p3Dopts.RotAngle = angleEarth(j);
planet3D('Earth Cloudy', p3Dopts);
hold on
p3Dopts = rmfield(p3Dopts, 'RotAngle');
sunRelPos = normalize(rSun(j,:),'norm')*dist';
sunRelSize = 2.1*Rsun/norm(rSun(j,:))*norm(sunRelPos-r);
[sunX,sunY,sunZ]=sphere;
sunX = sunX*sunRelSize+sunRelPos(1);
sunY = sunY*sunRelSize+sunRelPos(2);
sunZ = sunZ*sunRelSize+sunRelPos(3);
surf(sunX,sunY,sunZ,'EdgeColor','none','FaceColor',[0.98;0.843;0.627],'AmbientStrength',1,'SpecularStrength',1,'SpecularExponent',500)
%Sunlight, apparent Sun is ~28-34 arc minutes near Earth
light("Style","infinite","Position",sunRelPos)
light("Style","local","Position",sunRelPos+normalize(r-sunRelPos,'norm')*2*sunRelSize);
p3Dopts.Position = normalize(rMoon(j,:),'norm')*dist';
p3Dopts.Size = norm(p3Dopts.Position-r)/norm(rMoon(j,:)-r);
planet3D('Moon', p3Dopts);
% Flight path
scatter3( trackT(:,1), trackT(:,2), trackT(:,3), 12, scaledT, 'filled')
% Sun hitting the s/c
%scatter3( Y(:,1), Y(:,2), Y(:,3), 12, sunP, 'filled')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Earth equatorial frame', 'FontSize', 14);
cbar = colorbar; cbar.Color = 'w'; cbar.Title.Color = 'w';
cbar.Title.String = strcat('Time [',Tname,']');
clim([min(sunP);max(sunP)])   
axis equal;
grid on;
ax = gca; ax.Color = 'k'; ax.GridColor = 'w';
ax.GridAlpha = 0.45; ax.XColor = 'w'; ax.YColor = 'w';
ax.ZColor = 'w';
ax.CameraPosition = r;  % Set the camera positiont
%rotate3d;
ax.XLim = [-dist, dist];  % Set the x-axis range
ax.YLim = [-dist, dist];  % Set the y-axis range
ax.ZLim = [-dist, dist];  % Set the z-axis range
hold off
%camva(2*rad2deg(atan((norm(r)-Rearth)/norm(r))));

%% Functions
function r = relativeSun(tseconds, initialMJD, muSun)
T = initialMJD + tseconds/24/3600;
[kep,~] = uplanet(T, 3);
[r,~] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
r = -r;
end

function r = relativeMoon(tseconds, initialMJD)
T = initialMJD + tseconds/24/3600;
[r,~] = ephMoon(T);
r = r'; % 3x1
end