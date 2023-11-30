% Trajectory dynamics
function [ dY ] = dynamics(t, Y, propulsion, constants)

env = constants.env;
rocket = constants.rocket;

x = Y(1); vx = Y(2); z = Y(3); vz = Y(4); m = Y(5); gamma = Y(6);

% Changes in gravitational acceleration
g = env.g0*(env.R_E/(env.R_E+z))^2;

% Total velocity
V = sqrt(vx^2 + vz^2);

% Atmosphere
[rho,a,~,~,~,~,~] = atmos(z);
rho = 1.225 * exp(-z / 7.5e3);

% Compute mach
[~, a] = atmosisa(z);
M = V/a;

if t < 29.5
    CD = zerolift_drag_powered(M); %Positive thrust 
elseif vz > 0
    CD = zerolift_drag_coast(M) ;  %Ballistic trajectory
else
    if m<24
        CD = zerolift_drag_capsule(M);  %Rentry of the capsule + parachute TODO
    else 
        CD = 0.3;                       %TO DO Rentry of the launcher without paylaod
    end
end

% Calculate drag force
D = 0.5 * rho * CD * rocket.A * V^2;

% Thrust at time t
T = propulsion;
    
% Mass flow rate
m_dot=T/(rocket.Isp*env.g0);

% Gamma
gamma_dot = (-g/V+V/(env.R_E+z))*cos(gamma); 

% Calculate resultant force
Fx  =  T*cos(gamma) - D*cos(gamma);
Fz  =  T*sin(gamma) - D*sin(gamma) - m*g;

% Equations of motion with given thrust steering angle (theta)
dY = [
    vx;...
    Fx /m;...
    vz;...
    Fz / m;...
    -m_dot;...
    gamma_dot];
end