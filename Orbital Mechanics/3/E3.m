clear
close all

%% Input data
t1Window = [date2mjd2000([2023;11;1;0;0;0]); date2mjd2000([2025;1;1;0;0;0])];
t2Window = [date2mjd2000([2024;4;1;0;0;0]); date2mjd2000([2025;3;1;0;0;0])];
window1 = (t1Window(2)-t1Window(1))*24*3600;
window2 = (t2Window(2)-t2Window(1))*24*3600;

muSun = astroConstants(4);

p1 = 3;
p2 = 1;
porkchopUL = 50;

%% Porkchop plot and optimization
ngrid = 200;
[T1,T2] = meshgrid(linspace(t1Window(1),t1Window(2),ngrid),...
                   linspace(t2Window(1),t2Window(2),ngrid));
vCostgrid = ones(ngrid,ngrid)*porkchopUL;
deltaTgrid = zeros(ngrid,ngrid);
for i=1:ngrid
    for j=1:ngrid
        deltaTgrid(i,j) = T2(i,j)-T1(i,j);

        [vCost,vc1,vc2,eflag] = transferVcost(T1(i,j),T2(i,j),p1,p2,muSun);
        if (eflag == 0 && vc1<=7)
            vCostgrid(i,j) = min(vCost, porkchopUL);
        end
    end
end
% Best grid point
[Mp,Ip] = min(vCostgrid);
[M,I] = min(Mp);
point = [T1(Ip(I),I); T2(Ip(I),I)];

% Optimization
% Define function to solve for
func = @(x) transferVcost(x(1),x(2),p1,p2,muSun);

% Define constraint (times inside the departure and arrival spans)
lb = [t1Window(1); t2Window(1)];
ub = [t1Window(2); t2Window(2)];

% Solve
options = optimoptions('fmincon', 'Display', 'off');
[point, minVcost] = fmincon(func, point, [], [], [], [], lb, ub, [], options);
t1 = point(1); t2 = point(2);

deltaT = (t2-t1)*24*3600; %deltaT or time of flight (ToF) in [s]
[scaledT, Tname] = timescaling(deltaT);

% Porkchop
for i=1:ngrid
    for j=1:ngrid
        date = round(mjd20002date(T1(i,j)));
        T1(i,j) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
        date = round(mjd20002date(T2(i,j)));
        T2(i,j) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
    end
end
date = round(mjd20002date(t1));
point(1) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
date = round(mjd20002date(t2));
point(2) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));

figure;
% Contour of lines with equal times of flight
contour(T1,T2,deltaTgrid,'ShowText','on', 'LevelStep',60,'EdgeColor','black')
hold on
% Contour of curves with equal deltaV
contour(T1,T2,vCostgrid,'LineWidth',2)
% Lowest deltaV in the grid
scatter(point(1),point(2),30,'filled','red')
title(strcat('Lowest deltaV:', {' '}, num2str(minVcost), '[km/s]', {' '},...
    '(', num2str(scaledT), {' '}, Tname, ')'))

xlabel('Departure date'); ylabel('Arrival date');
cbar = colorbar;
clim([min(vCostgrid,[],'all'), porkchopUL]);
cbar.Title.String = strcat('deltaV [km/s]');
datetick('x',1,'keeplimits')
datetick('y',1,'keeplimits')
grid on;
hold off

rM = 0; % Prograde orbit. 1 if retrograde
Nrev = 0; % Number of revolutions
Ncase = 0; % Only when Nrev > 0. Small semi-major axis, or 1 for the large one

%% Get ephemerides
[kep, ~] = uplanet(t1Window(1), p1); %Earth
[r1Window, v1Window] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
[kep, ~] = uplanet(t1, p1); %Earth
[r1, v1] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
[kep, ~] = uplanet(t2, p1); %Earth
[r1p, v1p] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
[kep, ~] = uplanet(t2Window(1), p2); %Mars
[r2Window, v2Window] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
[kep, ~] = uplanet(t2, p2); %Mars
[r2, v2] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
[kep, ~] = uplanet(t1, p2); %Mars
[r2p, v2p] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');


%% Lambert problem
[a,p,e,eflag,vt1,vt2,deltaTparabolic,deltaTheta] = lambertMR( r1, r2, deltaT, muSun, rM, Nrev, Ncase );
vt1 = vt1'; vt2 = vt2'; %Velocity at points 1 and 2 of the transfer orbit

cost = vecnorm(vt1-v1) + vecnorm(v2-vt2);

%% Orbit propagation
% Initial state vector
y0 = [r1; vt1];

% Time grid
T = linspace( 0, deltaT, 1000 )';

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,muSun, 0, 0), T, y0, options );

% Scale time
[scaledT, Tname] = timescaling(T);

%% Plots
% Plot the flight path
figure;
% Low precision orbits: Initial
Yplot = timed2BP(r1Window, v1Window, muSun, 100, window1); %Departure window
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410, 0.2],'LineWidth',15);
hold on
Yplot = timed2BP(r1, v1, muSun, 100); %Earth orbit
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0, 0.4470, 0.7410],'LineWidth',3);
Yplot = timed2BP(r1, v1, muSun, 100, deltaT); %Earth motion during transfer
h3 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
Yplot = timed2BP(r2Window, v2Window, muSun, 100, window2); %Arrival window
h4 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.9290, 0.6940, 0.1250, 0.3],'LineWidth',15);
Yplot = timed2BP(r2, v2, muSun, 100); %Mars orbit
h5 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
Yplot = timed2BP(r2, v2, muSun, 100, -deltaT); %Mars motion during transfer
h6 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
Yplot = timed2BP(r1, vt1, muSun, 100); %Transfer orbit
h7 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0.4660, 0.6740, 0.1880],'LineWidth',3);
% Planets
planet3dOptions.Units = 'km';
planet3dOptions.Size = 20;
planet3D('Sun', planet3dOptions);
planet3dOptions.Size = 2000;
planet3dOptions.Position = r1;
planet3D('Earth', planet3dOptions);
planet3dOptions.Position = r2p;
planet3D('Mars', planet3dOptions);
planet3dOptions.FaceAlpha = 0.35;
planet3dOptions.Position = r1p;
planet3D('Earth', planet3dOptions);
planet3dOptions.Position = r2;
planet3D('Mars', planet3dOptions);
% Flight path
scatter3(Y(:,1), Y(:,2), Y(:,3), 6, scaledT) %Transfer path
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Flight path', 'FontSize', 14);
cbar = colorbar;
cbar.Title.String = strcat('Time [',Tname,']');
clim([min(scaledT);max(scaledT)])
lgnd = legend([h1,h2,h3,h4,h5,h6,h7],'Departure Window','Earth orbit',...
    'Earth motion during transfer','Arrival window','Mars orbit','Mars motion during transfer',...
    'Transfer orbit','Location','northoutside','NumColumns',3);
lgnd.FontSize = 12;
legend boxoff
axis equal;
grid on;
hold off