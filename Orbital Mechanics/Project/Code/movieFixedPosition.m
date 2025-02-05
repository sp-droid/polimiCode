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
tWindow = [date2mjd2000([2021;11;1;0;0;0]);0];
nsteps = 6000;

% Minimum value for this case: 600 steps per full period
% Takes 4.9s real time seconds per image, so 5400 steps in a night (7.5h)
% -> 7-8 periods, but in order to make it smooth i'm using 3.5 periods

%% Inputs
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.1; % Sidereal day
muSun = astroConstants(4);
Rsun = astroConstants(3);
TTsun = 5778;
muMoon = astroConstants(20);
CS = load('egm96/egm96_to360.ascii', '-ascii');
nmax = 360;
c = 299792.458; % Speed of light in km/s

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
wEarth = 2*pi/Tearth;
[A,B] = legendreAB(nmax);

%% Orbit
% Perform the integration
opts.RelTol = 1e-12;
opts.AbsTol = 1e-13;
opts.perturbShow = true;
[ Y, T ] = timed2BP(y0,muEarth,opts,nsteps,[],3);

% Define time scale, time window and time in mjd
[scaledT, Tname] = timescaling(T);

tWindow(2) = tWindow(1)+T(end)/24/3600;
Tmjd = linspace(tWindow(1), tWindow(2), nsteps)';

%% Celestial bodies
rSun = zeros(nsteps,3);
rMoon = zeros(nsteps,3);
angleEarth = zeros(nsteps,1);
track = zeros(nsteps,3);
for j=1:nsteps
    rSun(j,:) = relativeSun(T(j), tWindow(1), muSun);
    angleEarth(j) = rad2deg(wrapTo2Pi(wEarth*T(j)));
    track(j,1:3) = (rotRz(deg2rad(angleEarth(j)))*normalize(Y(j,1:3)','norm'))'*Rearth;
    rMoon(j,:) = relativeMoon(T(j), tWindow(1));
end

%% Animation
screen = get(0, 'ScreenSize');
v = VideoWriter ('videos/movieFixedPosition.avi');
open(v);
dist = max(vecnorm(Y(:,1:3)'));
r  = normalize([-2.7425;7.5133;-0.7491]','norm')*dist; %Get it with ax.CameraPosition
robs = r;
tic
for j=1:nsteps
    j
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
    p3Dopts = rmfield(p3Dopts, 'Position');
    p3Dopts = rmfield(p3Dopts, 'Size');
    % Flight path
    scatter3(Y(1:j,1), Y(1:j,2), Y(1:j,3), 12, T(1:j), 'filled')
    % Ground track
    scatter3(trackT(1:j,1), trackT(1:j,2), trackT(1:j,3), 12, T(1:j), 'filled')
    title(datestr(datetime(mjd20002date(tWindow(1)+T(j)/86400))), 'FontSize', 20,'Color','w');  
    axis equal;
    grid on;
    ax = gca; ax.Color = 'k'; ax.GridColor = 'k';
    ax.GridAlpha = 0; ax.XColor = 'k'; ax.YColor = 'k';
    ax.ZColor = 'k';
    
    ax.CameraPosition = r;  % Set the camera position
    camva(95.9378);
    ax.XLim = [-dist, dist];  % Set the x-axis range
    ax.YLim = [-dist, dist];  % Set the y-axis range
    ax.ZLim = [-dist, dist];  % Set the z-axis range

    hold off
    drawnow;

    frame = getframe (gcf);
    writeVideo (v, frame);
    close
end
toc
close(v);

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