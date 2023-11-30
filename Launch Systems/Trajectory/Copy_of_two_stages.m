clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Constants
graphPoints = 150;

% Time setup
t0 = 0;
tf = 1000;

[T0,tBurnout1] = Thrust_ASAS_13(0);
tHistBurnout1 = linspace(0, tBurnout1, graphPoints);

[T02,tBurnout2] = Thrust_Star13(0);
tHistBurnout2 = linspace(tBurnout1+2,tBurnout1+tBurnout2, graphPoints);
% Variables
zTurn = 500;

% Environmental constants
env.rhoRef = 1.225;                 % Air density at sea level
env.hScale = 7.5e3;                 % Height scale
env.g0 = 9.81;                      % Earth standard grav. acc.
env.R_E = 6278e3;                   % Earth standard radius

% Rocket constants
rocket.m0 = 250;                    % Rocket mass at t0 [kg]
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
%Y = [x; vx; z; vz; m; gamma];
Y0 = [0; 0.1; 0; 1; rocket.m0+payload.m; deg2rad(45)];


%% Thrust

thrust1 = arrayfun(@(t) Thrust_ASAS_13(t), tHistBurnout1);
thrust2 = arrayfun(@(t) Thrust_Star13(t), tHistBurnout2);

figure;
hold on
plot(tHistBurnout1,thrust1,'Color','blue','LineWidth',1.5)
plot(tHistBurnout2,thrust2,'LineWidth',1.5)
title(latex('Thrust vs time'),'Interpreter','latex');
xlabel(latex('Time [s]'),'Interpreter','latex');
ylabel(latex('Thrust [N]'),'Interpreter','latex');
grid on;
set(gca,'fontsize', 16)
hold off

%% Integration of trajectory

disp(['Rocket initial mass: ', num2str(Y0(5)), ' kg']);
disp(['Payload mass: ', num2str(constants.payload.m), ' kg']);
disp('--------------------')
% Gravity turn
options = odeset('Events', @(t, Y) zTurnEvent(t, Y, zTurn));
[tHist,Y]=ode78(@(t,y) rocketDynamics(t, y, false, constants), [t0,tf], Y0, options);
disp(['Gravity turn initiated at t+', num2str(tHist(end)), ' s']);
disp(['Gravity turn initiated at altitude: ', num2str(Y(end,3)), ' m'])
disp('---------')
% Stage separation 1

maxTimeStep = 0.5; % Maximum time step (adjust as needed)
options = odeset('Events', @stageSeparation1,'MaxStep', maxTimeStep);
Y0 = Y(end,:); Y0(5) = Y0(5)-25; % AMODIIIIIIFFIIERR
[tHist2,Y2]=ode78(@(t,y) rocketDynamics(t, y, true, constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist2]; Y = [Y; Y2];
disp(['Payload deployment at t+', num2str(tHist(end)), ' s']);
disp(['Payload deployment at altitude: ', num2str(Y(end,3)*1e-3), ' km']);

% Stage separation 2

maxTimeStep = 1; % Maximum time step (adjust as needed)
options = odeset('Events', @stageSeparation2,'MaxStep', maxTimeStep);
Y0 = Y(end,:); Y0(5) = Y0(5);
[tHist2,Y3]=ode78(@(t,y) rocketDynamics(t, y, true, constants), [tHist(end),tf], Y0, options);
tHist = [tHist; tHist2]; Y = [Y; Y3];
disp(['Payload deployment at t+', num2str(tHist(end)), ' s']);
disp(['Payload deployment at altitude: ', num2str(Y(end,3)*1e-3), ' km']);

% Rocket impact
options = odeset('Events', @impactEvent,'MaxStep', maxTimeStep);
Y0 = Y(end,:);Y0(5) = Y0(5)-constants.payload.m;
[tHist4,Y4]=ode78(@(t,y) rocketDynamics(t, y, true, constants), [tHist(end),tf-100], Y0, options);
tHist = [tHist; tHist4]; Y = [Y; Y4];

% % Payload trajectory
options = odeset('Events', @impactEvent,'MaxStep', maxTimeStep);
Y0 = Y3(end,:); Y0(5) = constants.payload.m;
[tHistP3,YP3]=ode78(@(t,y) rocketDynamics(t, y, true, constants), [tHist2(end),tf], Y0, options);
tHistP = tHistP3; YP = YP3;

figure;
plot(Y(:,1)*1e-3,Y(:,3)*1e-3,'Color',[0 0.4470 0.7410],'LineWidth',1.5,'DisplayName','Rocket')
hold on;
plot(YP(:,1)*1e-3,YP(:,3)*1e-3,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5,'DisplayName','Payload')
title(latex('Full trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
ylabel(latex('Altitude [km]'),'Interpreter','latex');
ylim([0, 1.05*max(Y(:,3))*1e-3]);
legend;
grid on;
set(gca,'fontsize', 16)
hold off;
%%
vx=Y(:,2);
vz=Y(:,4);
V=sqrt(vx.^2+vz.^2);
density = arrayfun(@(h) atmos(h), Y(:,3));
velocity = sqrt( Y(:,2).^2+Y(:,4).^2 );
dynamicPressure = 0.5*density.*V.^2;

maxiter = 1300;
% Velocity
figure;
yyaxis left
plot(Y(1:maxiter,1)*1e-3,V(1:maxiter),'LineWidth',1.5)
ylabel(latex('Velocity [m/s]'),'Interpreter','latex');
yyaxis right
plot(Y(1:maxiter,1)*1e-3,Y(1:maxiter,3)*1e-3,'LineWidth',1.5)
ylabel(latex('Altitude [km]'),'Interpreter','latex');
title(latex('Velocity vs trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
grid on;
xlim([0,max(Y(1:maxiter,1))*1e-3]);
ylim([0, 1.05*max(Y(1:maxiter,3))*1e-3]);
set(gca,'fontsize', 12)

% Dynamic pressure
figure;
yyaxis left
plot(Y(1:maxiter,1)*1e-3,dynamicPressure(1:maxiter)*1e-3,'LineWidth',1.5)
ylabel(latex('Pressure [kPa]'),'Interpreter','latex');
yyaxis right
plot(Y(1:maxiter,1)*1e-3,Y(1:maxiter,3)*1e-3,'LineWidth',1.5)
ylabel(latex('Altitude [km]'),'Interpreter','latex');
title(latex('Dynamic pressure vs trajectory'),'Interpreter','latex');
xlabel(latex('Downrange distance [km]'),'Interpreter','latex');
grid on;
xlim([0,max(Y(1:maxiter,1))*1e-3]);
ylim([0, 1.05*max(Y(1:maxiter,3))*1e-3]);
set(gca,'fontsize', 12)

% for i = 2:(length(V)-100)
%     a(i)=(V(i)-V(i-1))/(tHist(i)-tHist(i-1));
% end
% figure(3)
% plot(tHist(1:length(a)),a)
%% Functions
% Stops the function at impact, when z crosses 0
function [value, isterminal, direction] = impactEvent(t, Y)
    value = Y(3);                   % Triggers when z is 0
    isterminal = 1;                 % Stops the integration
    direction = -1;                 % Change in value must be decreasing
end
% Stage separation 1
function [value, isterminal, direction] = stageSeparation1(t, Y)
    % Define the event function
    [T0,tBurnout1] = Thrust_ASAS_13(0);
    value = tBurnout1-t;                   % Triggers when tburn1 is completed
    isterminal = 1;                 % Does not stop the integration
    direction = -1;                 % Change in value must be decreasing
end

% Stage separation at apogee
function [value, isterminal, direction] = stageSeparation2(t, Y)
    % Define the event function
    value = Y(4);                   % Triggers when dz is 0
    isterminal = 1;                 % Does not stop the integration
    direction = -1;                 % Change in value must be decreasing
end
% Stops the function at gravity turn point
function [value, isterminal, direction] = zTurnEvent(t, Y, zTurn)
    value = Y(3)-zTurn;             % Triggers when z is 0
    isterminal = 1;                 % Stops the integration
    direction = 0;                  % Change in value must be decreasing
end

function [ dY ] = rocketDynamics(t, Y, zTurn, constants)
[T02,tBurnout2] = Thrust_Star13(0);
[T0,tBurnout1] = Thrust_ASAS_13(0);
env = constants.env;
rocket = constants.rocket;
payload = constants.payload;

x = Y(1); vx = Y(2); z = Y(3); vz = Y(4); m = Y(5); gamma = Y(6);

g = env.g0*(env.R_E/(env.R_E+z))^2; %Gravity turn unit in [m]

% Total velocity
V = sqrt(vx^2 + vz^2);

% Compute mach 
M = computeMach(vx, vz, z);

if t < tBurnout2
    CD = zerolift_drag_powered(M); %Positive thrust 
elseif vz > 0
    CD = zerolift_drag_coast(M) ;  %Ballisitc trajectory
elseif vz < 0
    if m<70
        CD = zerolift_drag_capsule(M);  %Rentry of the capsule + parachute TODO
    else 
        CD = 0.4;                       %TO DO Rentry of the launcher without paylaod
    end
else 
    CD=0.3;
end



% Calculate air density
rho = env.rhoRef * exp(-z / env.hScale);

% Calculate drag force
Dx = 0.5 * rho * CD/2 * rocket.A *vx^2*cos(gamma);
Dz = 0.5 * rho * CD * rocket.A *vz^2*sin(gamma);


% Thrust at time t 
if t < tBurnout1
    T=Thrust_ASAS_13(t);
elseif t < tBurnout1+tBurnout2
    T=Thrust_Star13(t-tBurnout1);
else
    T=0;
end
    
%Mass flow rate
m_dot=T/(rocket.Isp*env.g0);

if zTurn==false
    gamma_dot = 0;
    Fx = 0;
    vx = 0;
    Fz = T-Dz-g*m;
else
    gamma_dot = -1/V*(g-V^2/(env.R_E+z))*cos(gamma); 
    
    % Calculate resultant force
    Fx  =  T*cos(gamma) - Dx;%*cos(gamma)
    Fz  =  T*sin(gamma) - Dz - m*g;
end


% Equations of motion with given thrust steering angle (theta)
dY = [
    vx;...
    Fx /m;...
    vz;...
    Fz / m;...
    -m_dot;...
    -gamma_dot];

end

function Mach = computeMach(vx, vz, z)
    [armand, vson] = atmosisa(z);
    Mach = sqrt(vx^2 + vz^2) / vson;
end


%% Miscellaneous functions

% Strings to latex
function string = latex(string)
    string = ['\boldmath$', string, '$'];
    string = strrep(string, ' ', '\hspace{0.5em}');
end





