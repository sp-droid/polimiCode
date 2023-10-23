clear all
close all
clc
 
%% DATA

% dry mass
m_dry = 400; % [kg]

% Total delta_V for maneuvers
DV_tot = 110; % [m/s]

% real margin - i.e. not the 100% that was taken during lecture
DV_margin = 91; % [m/s]
DV_tot = DV_tot + DV_margin; % [m/s]

% specific impulse
I_sp = (285 - 200) / (25 - 0.25) * 4.4 + 200;

% propellant density
rho_prop = 1.01 * 1e3; % [kg/m^3]

% Blowdown initial pressure [Pa]
P_i = 28.9 * 1e5; 

% Blowdown final pressure [Pa]
P_f = 5.2 * 1e5; 

% Helium specific gas constant
R_he = 2077.3; % [J / (kg*K)]

% Nitrogen specific gas constant
R_n2 = 296.80; % [J / (kg*K)]

% Initial Tank temperature
T_tank_i = 25.3 + 293.15; % [K]

% Titanium density and sigma
rho_Ti = 2780; % [kg / m^3]
sigma_Ti = 950 * 1e6; % [Pa]

% Earth gravity acceleration
g_0 = 9.81; % [m/s^2]

%% PROPELLANT

% final mass with 20% of margin
m_fin = m_dry;                                                        

% initial to final mass ratio
MR = exp(DV_tot / (g_0 * I_sp));

% initial mass 
m_in = MR*m_fin; % [kg]

% propellant mass with 5.5% margin
m_prop = m_in - m_fin; % [kg]
m_prop = m_prop * 1.055; 

% propellant volume with 10% margin
V_prop = m_prop / rho_prop; % [m^3]
V_prop = 1.10 * V_prop; 

%% HELIUM 
   
% Blow down ratio
B = P_i / P_f; 

% Initial volume of helium
V_in_gas = V_prop / (B - 1); % [m^3]

% Initial mass of helium 
m_he = (P_i * V_in_gas) / (R_he * T_tank_i) * 1.2; % [kg]


%% TANK

% tank volume with 1% of margin due to diaphgram
V_tank = 1.01 * (V_prop + V_in_gas); % [m^3]

% radius, thickness and mass of the tank
r_tank = ( 0.75 * V_tank / pi ) ^ (1/3); % [m]
t_tank = P_i * r_tank / (2 * sigma_Ti); % [m]
m_tank = rho_Ti * 4/3 * pi * ( (r_tank + t_tank)^3 - r_tank^3 ); % [kg]

%% RESULTS

RESULTS = [m_prop; m_he; m_dry; m_in; m_tank; r_tank; V_tank];

%% CHECK
% V_real = 0.12;
% r_real = ( 0.75 * 0.12 / pi ) ^ (1/3);
% m_real_prop = 22.3 + 17.5;
% m_real_he = 0.16;
% 
% err_r = abs(r_real - r_tank) / r_real;
% err_V = abs(V_tank - 0.12) / 0.12;
% err_m_prop = abs(m_real_prop - m_prop) / m_real_prop;
% err_m_he = abs(m_real_he - m_he) / m_real_he;




