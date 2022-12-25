clear
close all

%% Input data
mu_E = astroConstants(13);
rM = 0; % Prograde orbit. 1 if retrograde
Nrev = 0; % Number of revolutions
Ncase = 0; % Only when Nrev > 0. Small semi-major axis, or 1 for the large one
[r1,v1] = kep2car(12500,0,0,0,0,120,mu_E,'deg');
[r2,v2] = kep2car(9500,0.3,0,0,0,250,mu_E,'deg');
deltaT = 3300; %deltaT or time of flight (ToF) in [s]

%% Lambert problem
[a,p,e,eflag,vt1,vt2,deltaTparabolic,deltaTheta] = lambertMR( r1, r2, deltaT, mu_E, rM, Nrev, Ncase );
vt1 = vt1'; vt2 = vt2'; %Velocity at points 1 and 2 of the transfer orbit

cost = vecnorm(vt1-v1) + vecnorm(v2-vt2)

%% Orbit propagation
% Initial state vector
y0 = [r1; vt1];

% Perform the integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[ Y, T ] = timed2BP(y0,mu_E,opts,1000,deltaT);

% Scale time
[scaledT, Tname] = timescaling(T);

%% Plots
% Plot the flight path
figure;
% Black background
background('Black');
hold on
% Low precision orbits: Initial
Yplot = timed2BP([r1;v1], mu_E, [], 100);
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3)
Yplot = timed2BP([r2;v2], mu_E, [], 100);
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3)
Yplot = timed2BP([r1;vt1], mu_E, [], 100);
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',3)
hold on
% Planets
planet3dOptions.FaceAlpha = 1;
planet3dOptions.Units = 'km';
planet3D('Earth', planet3dOptions);
hold on
% Flight path
scatter3(Y(:,1), Y(:,2), Y(:,3), 6, scaledT)
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