clear all 
clc 
close all 

% Telecomunication subsystem 

k = 1.38*10^(-23); % [WS/k] Boltzmann constant 

% dara rates at Pluto estimation 
T_window = 172*8*3600; %[s]
Data = 5*10^9; %[Gbit]
R_new = Data / T_window; %[bps]


% Downlink = telemetry + scientific data 
f_downlink = 8.438*10^9; %[Hz] downlink frequency
%f_downlink = 12*10^9;
f_uplink = 7.2; %[GH] uplink frequency 
AU = 1.496*10^11; %[m]R_pluto = 4*10^3; %[bps] data rate at Pluto ( both telemetry and payload datas are included)--> the max is 4 kps
r_pluto = 36*AU; % average distance of Pluto wrt Sun
r_jupiter = 5.2*AU; % average distance of Jupiter wrt Sun
r_KBO = 50*AU; 
R_jupiter = 38*10^3; %[kbps] data rate at Jupiter
R_pluto = 2/6*10^3; %[bps] data rate at Pluto ( both telemetry and payload datas are included)--> the max is 4 kps
R_KBO = 5*10^9/(4.5*365*12*24*3600); %[bps] data rate at KBO 

time_Jupiter = 1; % time it took to reach Jupiter
time_Pluto = 2015-2006; % time it took to reach Pluto
time_KBO = 2016-2006; % time it took to reach KB0
Power_launch = 240; % [W] power at launch (RTG)
Power_decay = 3.5; % [W] over a period of 1 year; 
Power_Pluto = 200; % [W]
Power_Jupiter = Power_launch - time_Jupiter*Power_decay ;
Power_Pluto = Power_launch - time_Pluto*Power_decay;
Power_KB0 = Power_launch - time_KBO*Power_decay;


BER_downlink = 1e-5; 
BER_uplink = 1e-7;

Power_telecom_margin = 0.23; %[-] allocation of power to TTCS

% Antennas 
D_HGA = 2.1; %[m] diameter of the HGA 

% G_HGA = 42; %[decibels] minimum HGA gain 
mu_par = 0.55; %[-] efficiency of the antenna 
D_MGA = 0.3683; %[M] diameter of MGA 

f_downlink = 8.438*10^9; %[Hz] downlink frequency 
c = 3*10^8; %[m/s] speed of light
lambda_HGA = c/f_downlink; %[m] wavelength
% sizing og the antenna --> Gain 
A_HGA = pi*D_HGA^2 / 4; 
G_par_HGA = mu_par*4*pi*A_HGA / (lambda_HGA)^2; 
G_HGA = 10*log10(G_par_HGA);  %[db]
G_tx_HGA = G_HGA; 
theta_HGA = 65.3*lambda_HGA/ D_HGA;  % Beamwidth
margin = 1 - 10/100; %[deg] 
%theta_HGA = margin*theta_HGA; 
lambda_MGA = c/f_downlink; %[m] wavelength
% sizing og the antenna --> Gain 
A_MGA = pi*D_MGA^2 / 4; 
G_par_MGA = mu_par*4*pi*A_MGA / (lambda_MGA)^2; 
G_MGA = 10*log10(G_par_MGA);  %[db]
G_tx_MGA = G_MGA; 
theta_MGA = 65.3*lambda_MGA/ D_MGA;  % Beamwidth
margin = 1 - 10/100; %[deg] 
mu_GS = 0.01; 


% DS = uplink (it's receiving datas) 
T_GS = 21; % [K]
D_GS = 70 ; %[m] diameter of the DS antennas 
f_GS_downlink =8.420*10^9; %[Hz]
lambda_GS = c/f_GS_downlink; %[m] wavelength of the ground station
P_GS = 20000; %[W]
real_G_GS = 74.55; %[db] 
A_GS = pi*D_GS^2 / 4; % [m^2]
G_GS = mu_par*4*pi*A_GS / lambda_GS^2; 
G_GS = 10*log10(G_GS); % [db]
G_rx = real_G_GS; 

theta_GS= 65.3*lambda_GS/ D_GS  % Beamwidth
margin = 1; %[-] WORST case senario = 0.10 is lost 
theta_GS = margin*theta_GS;


% modulation = BPSK --> Salomon 
alpha_mod = 1; 
% CCSDs TURBO CODINDG
% it expanded by a factor of 6 the bandwidth 
alpha_enc = 6; 
delta_EB_turboenc = 3; %[db]
R_real = R_pluto*alpha_enc/alpha_mod; 
% thanks to turbo encoding 
EB_N0.min = 5.5- delta_EB_turboenc; %[dB] minimum Eb/N0 


P_tx= 24; 
P_input = P_tx / 0.53
P_tx = 10*log10(P_tx)

L_cable = -1 ; 
EIRP = P_tx + G_tx_HGA + L_cable
L_atm = -5*10^(-2);
L_space_pluto = 20*log10(lambda_HGA/(4*pi*r_pluto)) 
theta_GS = 65.3*lambda_GS/D_GS; 
L_point = -12*(mu_GS/theta_GS)^2; 
P_rx = EIRP+L_space_pluto+G_rx+L_point
R_data = 10*log10(R_real ) ; 
N0 = 10*log10(k*T_GS);
sizing.EB_N0_pluto = P_rx - R_data - N0 

% Sizing of carrier modulation index reduction 
beta_mod = deg2rad(30); %[deg]
P_mod_loss = 20*log10(cos(beta_mod)) %[db]
P_carrier = P_rx + P_mod_loss

% Signal to noise ratio 
SNR_carrier = P_carrier - N0 - 10*log10(10)
SNR_min = 10; %[db]
sizing.SNR_margin_pluto = SNR_carrier - SNR_min

Data = 5*10^9; %[bit]
T_windows = Data/(2*10^3*3600*24); 

