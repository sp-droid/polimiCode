clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Constants
graphPoints = 150;

% Time setup
t0 = 0;
tf = 1000;

[~,tBurnout1] = Thrust_ASAS_13(0);
tHistBurnout1 = linspace(0, tBurnout1, graphPoints);

[~,tBurnout2] = Thrust_Star13(0);
tHistBurnout2 = linspace(0,tBurnout2, graphPoints);

% Variables
zTurn = 500;

% Environmental constants
env.g0 = 9.81;                      % Earth standard grav. acc.
env.R_E = 6278e3;                   % Earth standard radius

% Rocket constants
rocket.m0 = 250;                    % Rocket EMPTY mass at t0 [kg]
rocket.CD = 0.4;                    % Rocket cd
rocket.D = 55.8*1e-2;               % Rocket diameter []
rocket.A = pi/4*rocket.D^2;         % Rocket area []
rocket.Isp = 283;                   % Rocket specific impulse [s]

% Payload constants
payload.m = 23;                     % Payload mass [kg]

constants.env = env;
constants.rocket = rocket;
constants.payload = payload;

% Initial conditions
%Y = [x; vx; z; vz; m];
Y0 = [0; 0; 0; 0; rocket.m0+payload.m];


%% Thrust

thrust1 = arrayfun(@(t) Thrust_ASAS_13(t), tHistBurnout1);
thrust2 = arrayfun(@(t) Thrust_Star13(t), tHistBurnout2);

figure;
hold on
plot(tHistBurnout1,thrust1*1e-3,'LineWidth',1.5,'DisplayName','ASAS 13')
plot(tHistBurnout2+tBurnout1,thrust2*1e-3,'LineWidth',1.5,'DisplayName','Star 15')
title(latex('Thrust vs time'),'Interpreter','latex');
xlabel(latex('Time [s]'),'Interpreter','latex');
ylabel(latex('Thrust [kN]'),'Interpreter','latex');
legend;
grid on;
set(gca,'fontsize', 12)

%% Integration of trajectory

disp(['Rocket initial mass: ', num2str(Y0(5)), ' kg']);
disp(['Payload mass: ', num2str(constants.payload.m), ' kg']);
disp('--------------------')
% From Launch to Gravity turn
options = odeset('Events', @(t, Y) zTurnEvent(t, Y, zTurn));
[tHist1,Y1] = ode78(@(t,y) rocketDynamics(t, y, Thrust_ASAS_13(t), constants), [t0,tf], Y0, options);
tHist = tHist1; Y = Y1;
disp(['Gravity turn initiated at t+', num2str(tHist(end)), ' s']);
disp(['Gravity turn initiated at altitude: ', num2str(Y(end,3)), ' m']);
disp('---------')
% From Gravity turn to 1st Burnout (assuming the turn happens before the first burnout)
options = odeset('Events', @(t, Y) burnoutEvent(t, Y, 0, tBurnout1));
Y0 = Y(end,:);
[tHist2,Y2] = ode78(@(t,y) rocketDynamics(t, y, Thrust_ASAS_13(t), constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist2]; Y = [Y; Y2];
disp(['First stage burnout at t+', num2str(tHist(end)), ' s']);
disp(['First stage burnout at altitude: ', num2str(Y(end,3)*1e-3), ' km']);
disp('---------')
% From 1st Burnout to 2nd Burnout (assuming the turn happens before the first burnout)
options = odeset('Events', @(t, Y) burnoutEvent(t, Y, tBurnout1, tBurnout2));
Y0 = Y(end,:);
[tHist3,Y3] = ode78(@(t,y) rocketDynamics(t, y, Thrust_Star13(t-tHist(end)), constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist3]; Y = [Y; Y3];
disp(['Second stage burnout at t+', num2str(tHist(end)), ' s']);
disp(['Second phase initiated at altitude: ', num2str(Y(end,3)*1e-3), ' km']);
disp('---------')
% From 2nd Burnout to Stage separation
options = odeset('Events', @stageSeparationEvent);
Y0 = Y(end,:);
[tHist4,Y4] = ode78(@(t,y) rocketDynamics(t, y, 0, constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist4]; Y = [Y; Y4];
disp(['Payload deployment at t+', num2str(tHist(end)), ' s']);
disp(['Payload deployment at altitude: ', num2str(Y(end,3)*1e-3), ' km']);
% Payload trajectory
% From Stage separation to Landing
options = odeset('Events', @impactEvent);
Y0 = Y(end,:); Y0(5) = constants.payload.m;
[tHist5,Y5]=ode78(@(t,y) rocketDynamics(t, y, 0, constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist5]; Y = [Y; Y5];

density = arrayfun(@(h) atmos(h), Y(:,3));
velocity = sqrt( Y(:,2).^2+Y(:,4).^2 );
dynamicPressure = 0.5*density.*velocity.^2;
acceleration = gradient(velocity);%sqrt( gradient(Y(:,2)).^2+gradient(Y(:,4)).^2 );
impactGs = velocity(end)/(tHist(end)-tHist(end-1))/constants.env.g0;

% Rocket trajectory after deploying payload
% From Stage separation to Impact
options = odeset('Events', @impactEvent);
Y0 = Y4(end,:); Y0(5) = Y0(5)-constants.payload.m;
[tHistR,YR]=ode78(@(t,y) rocketDynamics(t, y, 0, constants), [tHist4(end),tf], Y0, options);

velocityR = sqrt(YR(:,2).^2+YR(:,4).^2);
accelerationR = gradient(velocityR);%sqrt( gradient(YR(:,2)).^2+gradient(YR(:,4)).^2 );

disp('--------------------')
disp(['Burnt fuel: ', num2str(Y(1,5)-Y(end,5)), ' kg']);
disp(['Apogee: ', num2str(max(Y(:,3))*1e-3), ' km']);
disp(['Range: ', num2str(Y(end,1)*1e-3), ' km']);
disp(['Max speed: ', num2str(max(velocity)), ' m/s']);
disp(['Max acceleration: ', num2str(max(acceleration)/constants.env.g0), ' Gs']);
disp('--------------------')
disp(['Booster ground impact at t+', num2str(tHistR(end)), ' s']);
disp(['Impact velocity: ', num2str(velocityR(end)), ' m/s']);
disp(['Impact deceleration: ', num2str(velocityR(end)/(tHistR(end)-tHistR(end-1))/constants.env.g0), ' Gs']);
disp(['Impact energy: ', num2str(0.5*YR(end,5)*velocityR(end)^2/4.184*1e-6), ' kg of TNT']);

% Trajectory figure
figure;
plot(YR(:,1)*1e-3,YR(:,3)*1e-3,'LineWidth',1.5,'DisplayName','Rocket debris')
hold on;
plot(Y1(:,1)*1e-3,Y1(:,3)*1e-3,'LineWidth',1.5,'DisplayName','Launch->Turn')
plot(Y2(:,1)*1e-3,Y2(:,3)*1e-3,'LineWidth',1.5,'DisplayName','First stage')
plot(Y3(:,1)*1e-3,Y3(:,3)*1e-3,'LineWidth',1.5,'DisplayName','Second stage')
plot(Y4(:,1)*1e-3,Y4(:,3)*1e-3,'LineWidth',1.5,'DisplayName','Coast phase')
plot(Y5(:,1)*1e-3,Y5(:,3)*1e-3,'LineWidth',1.5,'DisplayName','Payload deployment->Landing')
title(latex('Full trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
ylabel(latex('Altitude [km]'),'Interpreter','latex');
xlim([0,max(YR(:,1))*1e-3]);
ylim([0, 1.05*max(Y(:,3))*1e-3]);
legend;
grid on;
set(gca,'fontsize', 12)
hold off;

% Velocity
figure;
yyaxis left
plot(Y(:,1)*1e-3,velocity,'LineWidth',1.5)
ylabel(latex('Velocity [m/s]'),'Interpreter','latex');
yyaxis right
plot(Y(:,1)*1e-3,Y(:,3)*1e-3,'LineWidth',1.5)
ylabel(latex('Altitude [km]'),'Interpreter','latex');
title(latex('Velocity vs trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
grid on;
xlim([0,max(Y(:,1))*1e-3]);
ylim([0, 1.05*max(Y(:,3))*1e-3]);
set(gca,'fontsize', 12)

% Acceleration
figure;
yyaxis left
plot(Y(:,1)*1e-3,acceleration/constants.env.g0,'LineWidth',1.5)
ylabel(latex('Acceleration [Gs]'),'Interpreter','latex');
yyaxis right
plot(Y(:,1)*1e-3,Y(:,3)*1e-3,'LineWidth',1.5)
ylabel(latex('Altitude [km]'),'Interpreter','latex');
title(latex('Acceleration vs trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
grid on;
xlim([0,max(Y(:,1))*1e-3]);
ylim([0, 1.05*max(Y(:,3))*1e-3]);
set(gca,'fontsize', 12)

% Dynamic pressure
figure;
yyaxis left
plot(Y(:,1)*1e-3,dynamicPressure*1e-3,'LineWidth',1.5)
ylabel(latex('Pressure [kPa]'),'Interpreter','latex');
yyaxis right
plot(Y(:,1)*1e-3,Y(:,3)*1e-3,'LineWidth',1.5)
ylabel(latex('Altitude [km]'),'Interpreter','latex');
title(latex('Dynamic pressure vs trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
grid on;
xlim([0,max(Y(:,1))*1e-3]);
ylim([0, 1.05*max(Y(:,3))*1e-3]);
set(gca,'fontsize', 12)

%% Functions

% Stops the function at gravity turn point
function [value, isterminal, direction] = zTurnEvent(t, Y, zTurn)
    value = Y(3)-zTurn;             % Triggers when z is 0
    isterminal = 1;                 % Stops the integration
    direction = 0;
end
% Stops the function at burnout
function [value, isterminal, direction] = burnoutEvent(t, Y, t0, tBurnout)
    value = t-t0-tBurnout;          % Triggers when z is 0
    isterminal = 1;                 % Stops the integration
    direction = 0;
end
% Stage separation at apogee
function [value, isterminal, direction] = stageSeparationEvent(t, Y)
    % Define the event function
    value = Y(4);                   % Triggers when dz is 0
    isterminal = 1;                 % Does not stop the integration
    direction = -1;                 % Change in value must be decreasing
end
% Stops the function at impact, when z crosses 0
function [value, isterminal, direction] = impactEvent(t, Y)
    value = Y(3);                   % Triggers when z is 0
    isterminal = 1;                 % Stops the integration
    direction = -1;                 % Change in value must be decreasing
end

% Trajectory dynamics
function [ dY ] = rocketDynamics(t, Y, propulsion, constants)

env = constants.env;
rocket = constants.rocket;
payload = constants.payload;

x = Y(1); vx = Y(2); z = Y(3); vz = Y(4); m = Y(5);

% Gravity turn unit in [m]
g = env.g0*(env.R_E/(env.R_E+z))^2;

% Total velocity
V = sqrt(vx^2 + vz^2);

% Atmosphere
[rho,a,~,~,~,~,~] = atmos(z);

% Compute mach
M = V/a;

% if t < 55
%     CD = zerolift_drag_powered(M); %Positive thrust 
% elseif vz > 0
%     CD = zerolift_drag_coast(M) ;  %Ballistic trajectory
% else
%     if m<70
%         CD = zerolift_drag_capsule(M);  %Rentry of the capsule + parachute TODO
%     else 
%         CD = 0.3;                       %TO DO Rentry of the launcher without paylaod
%     end
% end
CD = 0.4;

% Calculate drag force
D = 0.5 * rho * CD * rocket.A * V^2;

% Thrust at time t
T = propulsion;
    
% Mass flow rate
m_dot=T/(rocket.Isp*env.g0);

% Gamma gamma_dot = (-g/V+V/(env.R_E+z))*cos(gamma); 
if V==0
    gamma = pi/2;
else
    gamma = asin(vz/V);
end

% Calculate resultant force
Fx  =  T*cos(gamma) - D*cos(gamma);
Fz  =  T*sin(gamma) - D*sin(gamma) - m*g;

% Equations of motion with given thrust steering angle (theta)
dY = [
    vx;...
    Fx /m;...
    vz;...
    Fz / m;...
    -m_dot];
end

%% Miscellaneous functions

% Strings to latex
function string = latex(string)
    string = ['\boldmath$', string, '$'];
    string = strrep(string, ' ', '\hspace{0.5em}');
end