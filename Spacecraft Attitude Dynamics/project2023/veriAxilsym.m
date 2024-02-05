%% SPACECRAFT ATTITUDE DYNAMICS
clear;
clc;
close all;
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.3 0.4]);

%% PARAMS
% Inertia and initial angles
I = [0.06; 0.01; 0.01];
% omega_sc_0 = [0.45; 0.52; 0.55];
%I = [0.05; 0.05; 0.01];
omega_sc_0 = [0; 0; 0];
gamma_0 = [0;0;1];

I = diag(I);
Inv_I = pinv(I);

A_0 = eye(3);

rr_Esc_0 = [7070.6; 0; 0];
vv_Esc_0 = [0; 7.4338; 0];
mu_E = 398600;
R=6378.1363;                 %Earth's equatorial radius [km]
R_E_mean = 6371.01;              % Mean radius of Earth          [km]
R_E_eqtr = 6378.1363;            % Equatorial radius of Earth    [km]
J2 = 1.082626925638815e-3;       % Earth's second zonal coefficient [-]   
omega_E = 2*pi/(24*60*60);       % Earth's rotational velocity [rad/s]

Winter_solstice_UTC = [2023, 12, 22,  3, 27, 0];
theta_G_Winter_solstice = 5.615596868291755;
DAY = [2023,12,22,6,0,0];
t_elapsed = seconds(datetime(DAY)-datetime(Winter_solstice_UTC));
% RAAN of Greenwich meridian at DAY
theta_G0 = mod(theta_G_Winter_solstice+omega_E*t_elapsed,2*pi);    % [rad]

% On the col diff m, on the row diff n

g_nm = [-29404.8 -1450.9     0      0     0
         -2499.6  2982.0  1677.0    0     0
          1363.2 -2381.2  1236.2  525.7   0
           903.0   809.5    86.3 -309.4  48.0];

h_nm = [  0  4652.5      0      0       0
          0 -2991.6  -734.6     0       0 
          0   -82.1   241.9  -543.4     0
          0   281.9  -158.4   199.7  -349.7];

% Coefficients to adjust yearly the gaussians
correct_g = [5.7   7.4   0      0     0
             -11   -7    -2.1   0     0
              2.2  -5.9  3.1   -12    0
             -1.2  -1.6  -5.9   5.2   -5.1];

correct_h = [  0  -25.9     0       0      0
               0  -30.2    -22.4    0      0 
               0   6       -1.1     0.5    0
               0  -0.1      6.5     3.6   -5];

smith_norm=zeros(4,5);

for n=1:4
    for m=0:4
        if(m==0)
        smith_norm(n,m+1)=(factorial(n-m)/factorial(n+m))^(1/2);
        elseif(n>=m)
        smith_norm(n,m+1)=(2*factorial(n-m)/factorial(n+m))^(1/2);
        else
        smith_norm(n,m+1)=0;
        end
    end
end

M_rw = 0.2;                              % mass of the single rw [kg]
M_c = 5.5+3*M_rw;                        % mass of the central body [kg]
M_a = 0.150;                             % mass of a single antenna [kg]
M_p = 0.650;                             % mass of a single solar array [kg]
M_tot = M_c+6*M_p+4*M_a;                 % total mass of the spacecraft
m_sc_body = [1;1;1]*10e-3*M_tot/sqrt(3);                 % [A*m^2]

%% SIMULATION
disp('Running simulink model...');
% Params
simfile = 'verificationGG.slx';

% Execution
out = sim(simfile);
time = out.tout;
w = out.omega_sc;
gamma = out.gamma*180/pi;
% gamma = out.gamma;

%% Post processing
% Calculate h and T
I = diag(I);
T = 0.5*(I(1)*w(:,1).^2+I(2)*w(:,2).^2+I(3)*w(:,3).^2);
T = T - T(1);
h = sqrt((I(1)*w(:,1)).^2+(I(2)*w(:,2)).^2+(I(3)*w(:,3)).^2);
h = h - h(1);

disp('Plotting graphs...');

%% Plots
doublePlot(time,T,time,h,'Time [s]','Kinetic energy T-T0 [J]','Angular momentum H-H0 [kg m^2/s]','Torque free, spinning around stable axis');
% simplePlot(time,{gamma(:,1),gamma(:,2),gamma(:,3)},'Time [s]','Magnitude []','Pointing direction',{'gammax','gammay','gammaz'});
simplePlot(time,{w(:,1),w(:,2),w(:,3)},'Time [s]','Angular velocity [rad/s]','Torque free, spinning around stable axis',{'wx','wy','wz'});
simplePlot(time,{gamma(:,1),gamma(:,2),gamma(:,3)},'Time [s]','Angle [deg]','LVLH angles',{'yaw','pitch','roll'});


% Known case
lambda = (I(3)-I(1))*omega_sc_0(3)/I(1);
w(:,3) = ones(1,length(time))*omega_sc_0(3);
w(:,2) = omega_sc_0(1)*sin(lambda*time)+omega_sc_0(2)*cos(lambda*time);
w(:,1) = omega_sc_0(1)*cos(lambda*time)-omega_sc_0(2)*sin(lambda*time);
%simplePlot(time,{w(:,1),w(:,2),w(:,3)},'Time [s]','Angular velocity [rad/s]','Torque free, axisymmetric, analytical',{'wx','wy','wz'});