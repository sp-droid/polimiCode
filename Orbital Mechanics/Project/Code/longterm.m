clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Inputs
% Physical parameters
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
norbs = 1464;
ngrid = 14640000;

[a,e,i,bOmega,sOmega,theta] = car2kep(r0,v0,muEarth,'rad');
kep0 = [a,e,i,bOmega,sOmega,theta]';

%% Orbit propagation
densityModel = @(r) densitySimplified(norm(r)-Rearth);
opts.j2Pert = @(r) j2Pert(r,J2,Rearth,muEarth);
opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, 2.1, 0.0095);

opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[Y1, T] = timed2BP(y0, muEarth, opts, ngrid, [], norbs);

[scaledT, Tname] = timescaling(T);

%% Filtering
% We want to remove all orbits but one each year or period
Tsize = 365.25*24*3600/12;
orbSize = T(end)/norbs*1.15;
k = 0;
for j=1:ngrid
    if (T(j)-k*Tsize>Tsize)
        k = k+1;
    elseif (T(j)-k*Tsize>orbSize)
        Y1(j,:) = [NaN NaN NaN NaN NaN NaN];
    end
end

%% Plots
[dist,j] = max(vecnorm(Y1(:,1:3)'));
r  = normalize([1.6317e+05;-80262;5.9325e+05]','norm')*dist; %Get it with ax.CameraPosition

figure('Color','k');
p3Dopts.Units = 'km';
planet3D('Earth', p3Dopts);
hold on
scatter3( Y1(:,1), Y1(:,2), Y1(:,3), 6, scaledT, 'filled')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
cbar = colorbar; cbar.Color = 'w'; cbar.Title.Color = 'w';
cbar.Title.String = strcat('Time [',Tname,']');
clim([min(scaledT);max(scaledT)])
axis equal;
grid on;
ax = gca; ax.Color = 'k'; ax.GridColor = 'w';
ax.GridAlpha = 0.25; ax.XColor = 'w'; ax.YColor = 'w';
ax.ZColor = 'w';
hold off