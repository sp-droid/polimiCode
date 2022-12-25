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

% Perform the integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[ Y, T ] = timed2BP(y0,mu_E,opts,1000,deltaT);

% Scale time
[scaledT, Tname] = timescaling(T);

%% Plots
% Plot the flight path
figure;
% Black background and planets
background('Black');
hold on
planet3dOptions.FaceAlpha = 1;
planet3dOptions.Units = 'km';
planet3D('Earth', planet3dOptions);
hold on
% Flight path
scatter3( Y(:,1), Y(:,2), Y(:,3), 6, scaledT)
h1 = plot3(Y(1,1), Y(1,2), Y(1,3),'^','Color',[0,1,0],'LineWidth',6);
h2 = plot3(Y(end,1), Y(end,2), Y(end,3),'v','Color',[0,0.8,0],'LineWidth',6);
hold on
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Flight path', 'FontSize', 14);
cbar = colorbar;
cbar.Title.String = strcat('Time [',Tname,']');
clim([min(scaledT);max(scaledT)])
lgnd = legend([h1,h2],'Start','End','Location','southwest','NumColumns',2);
lgnd.FontSize = 14;
legend boxoff
axis equal;
grid on;
ax = gca;
ax.GridColor = [1,1,1];
ax.GridAlpha = 0.25;
hold off