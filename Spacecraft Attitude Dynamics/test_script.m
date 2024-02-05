clear all
close all
clc

%% DATA
% earth orbit

ra = 3e4; rp = 8e3;
a = (ra+rp)/2;
e = (ra-rp)/(ra+rp);

orbit_E = struct('mu', astroConstants(13), 'a', a, 'e', e, ...
                 'inc', deg2rad(10), 'G', astroConstants(1), ...
                 'theta0', 0, 'R', astroConstants(23));

orbit_E.n = sqrt(orbit_E.mu/orbit_E.a^3);
orbit_E.T = 2*pi / orbit_E.n;
orbit_E.w = 2*pi/(23*3600+56*60+4.09);   

% solar orbit
orbit_S = struct('mu', astroConstants(4), 'a', astroConstants(2), 'e', 0.0167, ...
                 'epsilon', deg2rad(23.45), 'G', astroConstants(1), ...
                 'theta0', 0  );

orbit_S.n =  sqrt(orbit_S.mu/orbit_S.a^3);


radiation = 1358; %(620 albedo@h +157 earth@h)
density = 2.72e-12;

m_t = 5.9722e24; % earth mass

% Spacecraft (inertia, initial angle, torque)
cubesat = struct('m', .8 * 12, 'a', .2, 'b', .2, 'c', .3405);

Ix = 1/12 * cubesat.m * (cubesat.b^2 + cubesat.c^2);
Iy = 1/12 * cubesat.m * (cubesat.a^2 + cubesat.c^2);
Iz = 1/12 * cubesat.m * (cubesat.a^2 + cubesat.b^2);

I = [Ix; Iy; Iz];
I_inv = 1./I;
w0 = orbit_E.n* [1; 1; 1];

% Reaction wheels
RW = struct('M_max', 0.23e-3);

% Solar Panels and SRP geometric parameters

sp_width  = .57;                % m
sp_length = .57;                % m
sp_a   = sp_length * sp_width;  % m^2
cs_a1 = cubesat.a * cubesat.c;  % m^2
cs_a2 = cubesat.a * cubesat.b;  % m^2

solarPanel = struct('P', radiation/(3e8), 'w', sp_width, 'l', sp_length, ...
                    'rho_d', .1, 'rho_s', [.5, .5, .5, .5, .5 .5, .1, .1, .1, .1], ...
                    'area', [cs_a1 cs_a1 cs_a2 cs_a1 cs_a1 cs_a2 sp_a sp_a sp_a sp_a]);

pos_CG = .15*[1;1;1];                               % Displacement of CG from geometric center
n1=[1;0;0]; n2=[0;1;0]; n3=[0;0;1];

pos(:,1) = cubesat.b/2 *n1; pos(:,4) = -pos(:,1);
pos(:,2) = cubesat.a/2 *n2; pos(:,5) = -pos(:,2);
pos(:,3) = cubesat.c/2 *n3; pos(:,6) = -pos(:,3);

pos(:,7) = (cubesat.a + solarPanel.l)/2 *n2; pos(:,8) = pos(:,7);
pos(:,9) = - pos(:,7); pos(:,10) = pos(:,9);

solarPanel.r = pos-pos_CG;

% Magnetic perturbation
m_par = [0.01;0.05;0.01];
mag_tilt = deg2rad(11.5);
H0 = sqrt((-29615e-9)^2+(-1728e-9)^2+(5186e-9)^2);


% Initial DCM
A0=[1 0 0; 0 1 0; 0 0 1];

% Sensor
gyro = struct('c', 15e-3, 'ARW', deg2rad(0.15), 'RRW', deg2rad(0.0003), ...
              'freq', 100, 'accuracy', 2.4e-9); % the accurancy is taken from wikipedia [rad]
gyro.A_eps = [1, 0.2e-3, 0.2e-3;
                0.2e-3, 1, 0.2e-3
                0.2e-3, 0.2e-3, 1];

% Attitude Determination
sun_sensor    = struct('dvst', 3*0.2 ,'accuracy', deg2rad(0.2));
earth_horizon = struct('dvst', 3*0.2 ,'accuracy', deg2rad(0.2));

% PD parameters
K_p = Iz * (20 * orbit_E.n)^2;   % we assumed w = 20*n
K_d = 40 * Iz * 0.7 * 20 * orbit_E.n; % we assumed csi = 0.3 

k1 = .09;
k2 = .0085;

 
M_RW = 60e-3;
R_RW = 23e-3;
Ir_RW = .25 * M_RW * R_RW^4;
Ir = Ir_RW;

time = struct('de_tumble', (orbit_E.T)/100, 'traj_track', (orbit_E.T)/2);

%% COMPUTATION
w    = out.omega;
M_RW = out.M_RW;
M_magn = out.M_magn;
M_srp = out.M_srp;
M_grav = out.M_grav;
point_err = out.point_err;

time = 0:length(w);
%% PLOT

figure()
hold on
plot(time, w(1));
plot(time, w(2));
plot(time, w(3));

xlabel('time [s]'); ylabel('angular velocity [rad/s]');
legend('w_x', 'w_y', 'w_z');
title('Spacecraft angular velocity in body axes');

grid on;


