clear
close all

%% Input data
rM = 0; % Prograde orbit. 1 if retrograde
Nrev = 0; % Number of revolutions
Ncase = 0; % Only when Nrev > 0. Small semi-major axis, or 1 for the large one
r1 = [-21800; 37900; 0]; %[km]
r2 = [27300; 27700; 0]; %[km]
deltaT = 15*3600+6*60+40; %deltaT or time of flight (ToF) in [s]

mu_E = astroConstants(13);

%% Lambert problem
[a,p,e,eflag,v1,v2,deltaTparabolic,deltaTheta] = lambertMR( r1, r2, deltaT, mu_E, rM, Nrev, Ncase );
v1 = v1'; v2 = v2';

%% Orbit propagation
% Initial state vector
y0 = [r1; v1];

% Time grid
T = linspace( 0, deltaT, 1000 )';

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E, 0, 0), T, y0, options );

% Scale time
[scaledT, Tname] = timescaling(T);

%% Plots
% Plot the orbit
figure;
background('Black');
hold on
planet3dOptions.Color = [1;1;1];
planet3dOptions.Units = 'km';
planet3D('Earth', planet3dOptions);
hold on
scatter3( Y(:,1), Y(:,2), Y(:,3), 6, scaledT)
hold on
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Flight path');
cbar = colorbar;
cbar.Title.String = strcat('Time [',Tname,']');
clim([min(scaledT);max(scaledT)])
axis equal;
grid on;
ax = gca;
ax.GridColor = [1,1,1];
ax.GridAlpha = 0.25;
hold off