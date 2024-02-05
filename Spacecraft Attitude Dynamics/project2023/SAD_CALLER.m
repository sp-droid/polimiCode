%% SPACECRAFT ATTITUDE DYNAMICS
clear;
clc;
close all;
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.3 0.4]);
%%

% First select maneuver, then run to end, then run the simulation on
% simulink

maneuver=input('type of maneuver, 1 for nadir pointing, 2 for de-tumbling, 3 for slew ')

executionType = input ('Type of execution: 0 for matlab code (to run simulink separately), 1 for matlab + simulink + plots')

%% Celestial parameters initialization

% Celestial bodies gravitational constant

mu_S = 1.327124400e11;           % Sun gravitational constant   [km^3/s^2]
mu_E = 3.986004418e5;            % Earth gravitational constant [km^3/s^2]

% Earth's geometry

R_E_mean = 6371.01;              % Mean radius of Earth          [km]
R_E_eqtr = 6378.1363;            % Equatorial radius of Earth    [km]
J2 = 1.082626925638815e-3;       % Earth's second zonal coefficient [-]   
omega_E = 2*pi/(24*60*60);       % Earth's rotational velocity [rad/s]

% Earth's magnetic field

% IGRF table for the gaussians at year 2020
% report by National Centers for Environmental Information, USA

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

% Smith normalization to normalize the polynomials

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

% Time reference

Winter_solstice_UTC = [2023, 12, 22,  3, 27, 0];
theta_G_Winter_solstice = 5.615596868291755;

Summer_solstice_UTC = [2024,  6, 20, 20, 51, 0];

% Select time from which propagate the orbit

%DAY=Summer_solstice_UTC;
DAY = [2023,12,22,6,0,0];

t_elapsed = seconds(datetime(DAY)-datetime(Winter_solstice_UTC));


%% Spacecraft properties 

% MASS


% 6U Cubesat mass allocation.
% The spacecraft is divided into its basic components:
% Casing and payload
% 3 Reaction Wheels
% 2 Solar panels, each composed of 3 arrays
% 4 SAR antennas

M_rw = 0.2;                              % mass of the single rw [kg]
M_c = 5.5+3*M_rw;                        % mass of the central body [kg]
M_a = 0.150;                             % mass of a single antenna [kg]
M_p = 0.650;                             % mass of a single solar array [kg]
M_tot = M_c+6*M_p+4*M_a;                 % total mass of the spacecraft




% GEOMETRIC PROPERTIES 


% Central body nominal dimensions

x_c = 0.3;                               % length along the x-body-axis of the central body [m]                             
z_c = 0.2;                               % length along the z-body-axis of the central body [m]   
y_c = 0.1;                               % length along the y-body-axis of the central body [m] 

% Solar panels nominal dimensions

x_p = 0.2;                               % length along the x-body-axis of the panel [m]
y_p = 0.3;                               % length along the y-body-axis of the panel [m]
z_p = 0.002;                             % length along the z-body-axis of the panel [m]




% CENTER OF MASS


% Position of the center of mass:
% The spacecraft is symmetrical along the x-body axis, so the CM is located
% along the axis. The poisition in computed from the bottom base.

d_cc = 0.15;            % geometric center along x of the central body [m]

x_cm = ((d_cc+0.1)*M_c+(x_c+0.1)*6*M_p+0.1*M_a)/M_tot-0.1; % [m]

% Distances along x from the center of mass of the geometrical centers of the area
d_cc_cm = (x_cm-d_cc);  % [m]




% INERTIA OF THE COMPONENTS


% Inertia of the solar panels, centered at their geometric center

I_x_p = 1/12*M_p*(y_p^2+z_p^2);
I_y_p = 1/12*M_p*(z_p^2+x_p^2);
I_z_p = 1/12*M_p*(y_p^2+x_p^2);
I_p = diag([I_x_p,I_y_p,I_z_p]);
d_p = 0.2;
d_p2 = 0.5;
d_p3 = 0.8;
I_p = 2*(I_p+diag([M_p*d_p^2,0,M_p*d_p^2])+I_p+diag([M_p*d_p2^2,0,M_p*d_p2^2])+I_p+diag([M_p*d_p3^2,0,M_p*d_p3^2]));    % inertia of the panels centred on the tip o the spacecraft

% Inertia of the central body, centered at its geometric center

I_x_c = 1/12*M_c*(z_c^2+y_c^2);
I_y_c = 1/12*M_c*(x_c^2+z_c^2);
I_z_c = 1/12*M_c*(x_c^2+y_c^2);
I_c = diag([I_x_c,I_y_c,I_z_c]);         % inertia of the central body [kg*m^2]

% Inertia of a single antenna, centered at its geometric center

r_a = 0.020;                             % radius of the antenna [m]
y_a = 1.5;                               % length of the antenna [m]
I_a_c = 1/12*M_a*y_a^2;                  % radial intertia of one antenna [kg*m^2]
I_a_t = 1/2*M_a*r_a^2;                   % axial intertia of one antenna [kg*m^2]
d_a = 0.8;                               % distance from the center of the bottom part of the spacecraft [m]

% Inertia of the 2 antennas directed along the y-body-axis

I_a_y = diag([I_a_c,I_a_t,I_a_c]);       % [kg*m^2]

% Inertia of the 2 antennas directed along the z-body-axis

I_a_z = diag([I_a_c,I_a_c,I_a_t]);       % [kg*m^2]

% Inertia of the antennas, centred at the bottom of the spacecraft 

I_a = 2*(I_a_y+diag([M_a*d_a^2,0,M_a*d_a^2])+I_a_z+diag([M_a*d_a^2,M_a*d_a^2,0])); % [kg*m^2]   




% FULL SPACECRAFT INERTIA MATRIX


I = I_a+diag([0,4*M_a*x_cm^2,4*M_a*x_cm^2])+I_c+diag([0,M_c*d_cc_cm^2,M_c*d_cc_cm^2])+I_p+diag([0,6*M_p*(d_cc+d_cc_cm)^2,6*M_p*(d_cc-d_cc_cm)^2]);    % total inertia of the spacecraft
Inv_I = inv(I);

K_yaw=(I(3,3)-I(2,2))/I(1,1);
K_roll=(I(3,3)-I(1,1))/I(2,2);
K_pitch=(I(2,2)-I(1,1))/I(3,3);




% SURFACE MODELLING 


% Body faces normal unit vectors

N_1 = [0;0;1];                      
N_2 = [0;1;0];
N_3 = [0;0;-1];
N_4 = [0;-1;0];
N_5 = [1;0;0];
N_6 = [-1;0;0];

% Solar panels normal unit vectors

N_7 = N_1;
N_8 = N_1;
N_9 = -N_1;
N_10 = -N_1;
N = [N_1,N_2,N_3,N_4,N_5,N_6,N_7,N_8,N_9,N_10];

% Areas of the spacecraft faces

AA = [0.03;0.06;0.03;0.06;0.02;0.02;0.18;0.18;0.18;0.18];     % [m^2]

% Body faces vector distances from center of mass to geometrical center of
% area

d_spc_1 = 0.1*N_1+d_cc_cm*[-1;0;0];
d_spc_2 = 0.05*N_2+d_cc_cm*[-1;0;0];
d_spc_3 = 0.1*N_3+d_cc_cm*[-1;0;0];
d_spc_4 = 0.05*N_4+d_cc_cm*[-1;0;0];
d_spc_5 = 0.15*N_5+d_cc_cm*[-1;0;0];
d_spc_6 = 0.15*N_6+d_cc_cm*[-1;0;0];

% Solar panels vector distances from center of mass to geometrical center of
% area

d_spc_7 = 0.06*N_5+0.5*N_2+d_cc_cm*[-1;0;0];
d_spc_8 = 0.06*N_5+0.5*N_4+d_cc_cm*[-1;0;0];
d_spc_9 = 0.06*N_5+0.5*N_2+d_cc_cm*[-1;0;0];
d_spc_10 = 0.06*N_5+0.5*N_4+d_cc_cm*[-1;0;0];

d_spc = [d_spc_1,d_spc_2,d_spc_3,d_spc_4,d_spc_5,d_spc_6,d_spc_7,d_spc_8,d_spc_9,d_spc_10];




% SPACECRAFT RESIDUAL MAGNETIC DIPOLE


% No inspection on cubesat residual dipole
% Non spinning spacecraft
% Guidelines from NASA SP-8018

m_sc_body = [1;1;1]*10e-3*M_tot/sqrt(3);                 % [A*m^2]


%% Reaction wheel inertia matrix


% Inertia properties of a single reaction wheel
      
I_r =9.55e-5;               % inertia along spin axis [kg*m^2]

% Mounting of the RW ensamble wrt body frame

t_1=0.584196501093179;
t_2=-0.385718098741334;
t_3=0.374479392778978;

R = [0.9631    0.2014   -0.1783
    -0.2580    0.8792   -0.4005
     0.0761    0.4317    0.8988];

R_x=R(:,1);
R_y=R(:,2);
R_z=R(:,3);
v=[sqrt(3)/3,sqrt(3)/3,sqrt(3)/3];
R_rad_x_1=cross(R_x,v)/vecnorm(cross(R_x,v));
R_rad_x_2=cross(R_x,R_rad_x_1)/vecnorm(cross(R_x,R_rad_x_1));
RR_x=[R_x,R_rad_x_1',R_rad_x_2'];
R_rad_y_1=cross(R_y,v)/vecnorm(cross(R_y,v));
R_rad_y_2=cross(R_y,R_rad_y_1)/vecnorm(cross(R_y,R_rad_y_1));
RR_y=[R_y,R_rad_y_1',R_rad_y_2'];
R_rad_z_1=cross(R_z,v)/vecnorm(cross(R_z,v));
R_rad_z_2=cross(R_z,R_rad_z_1)/vecnorm(cross(R_z,R_rad_z_1));
RR_z=[R_z,R_rad_z_1',R_rad_z_2'];

R_xy=R(:,[1,2]);
R_yz=R(:,[2,3]);
R_xz=R(:,[1,3]);

inv_R = inv(R);
pinv_R_xy=pinv(R_xy);
pinv_R_yz=pinv(R_yz);
pinv_R_xz=pinv(R_xz);

M_brake=0.0002;

d_cmsc_rw = 0.003;                              % [m]
d_cm_rw = 0.00000005;                           % [m]


%% Horizon sensor coefficients

FOV = 34;                             % Sensor head field of view [deg]
Res = 1;                              % Resolution [deg]
f = (FOV/2)/tan(deg2rad(FOV/2));      % Normalized focal length of the sensor [deg]
alpha_mount_hs = 0.15;                % Mounting error in each direction [deg]


%% Sun sensor coefficients

alpha_mount_ss = -0.5;                 % Mounting error in each direction [deg]


%% Sun orbit paramaters

% Orbital parameters of the Sun apparent orbit in ECI frame as of 
% Winter Solstice 2023

a_S = 149597870;                % Semi-major axis [km]
e_S = 0.0167133;                % Eccentricity [-]
i_S = deg2rad(23.4406);         % Inclination [rad]
OM_S = 0;                       % RAAN [rad]
om_S = deg2rad(282.7685);       % Arg. of Pericenter [rad]
theta_S_start = 6.0603;         % True anomaly as of Winter Solstice 2023

% Poisition and velocity vectors of the Sun in ECI frame as of DAY

[theta_S_0] = mod(orbit_propagation(a_S,e_S,theta_S_start,t_elapsed,mu_S),2*pi);
[rr_ES_0,vv_ES_0] = koe2rv(a_S,e_S,i_S,OM_S,om_S,theta_S_0,mu_S); % [km;km/s]


%% Earth parameters

% RAAN of Greenwich meridian at DAY
theta_G0 = mod(theta_G_Winter_solstice+omega_E*t_elapsed,2*pi);    % [rad]         


%% Spacecraft orbit parameters

% Orbital parameters of the spacecraft in ECI frame as of 
% Winter solstice 2023

h = 702;                      % Altitude [m]
a = R_E_mean+h;               % Semi-major axis [m]
e = 0;                        % Eccentricity [-]
i = deg2rad(98.2);            % Inclination [rad]
OM = deg2rad(0);              % RAAN [rad]
om = 0;                       % Arg. of Perigee [rad]
theta = 2.8187;               % True anomaly [rad]
T = 2*pi*sqrt(a^3/mu_E);      % Orbital period [s]
n = 2*pi/T;                   % Mean angular velocity [rad/s]

% Spacecraft Position and velocity vectors from orbital parameters as of 
% Winter Solstice

[rr_Esc_start,vv_Esc_start] = koe2rv(a,e,i,OM,om,theta,mu_E); % [km;km/s]

% Orbit propagation until DAY, include J2 perturbation

options = odeset(RelTol=1e-13,AbsTol=1e-14);

[~,s] = ode113(@twobody_J2,[0,t_elapsed],[rr_Esc_start;vv_Esc_start],options,mu_E);
rr_Esc_0 = s(end,(1:3))';                                      % [km]
vv_Esc_0 = s(end,(4:6))';                                      % [km/s]


%% Attitude initial condition

pitch_Esc_0 = cross(rr_Esc_0,vv_Esc_0);
roll_Esc_0 = cross(pitch_Esc_0,rr_Esc_0);

LVLH_X_YAW_0 = (rr_Esc_0./vecnorm(rr_Esc_0))';
LVLH_Y_ROLL_0 = (roll_Esc_0 ./vecnorm(roll_Esc_0 ))'; 
LVLH_Z_PITCH_0 = (pitch_Esc_0./vecnorm(pitch_Esc_0))';

A_LVLH_ECI_0 = [LVLH_X_YAW_0;LVLH_Y_ROLL_0;LVLH_Z_PITCH_0];

switch maneuver
    case 1
        omega_sc_0=[10e-6,10e-6,n+10e-6];
        A_B_LVLH=[1,0.009258,0.0078;-0.0092,0.999,-0.006;-0.007693,-0.00576,0.999];
        A_0=A_B_LVLH*A_LVLH_ECI_0;
        A_slew_B_LVLH=zeros(3);             % unused in this maneuver

    case 2
        omega_sc_0=[1e-2,1e-2,1e-2];
        A_B_LVLH=[1,0.009258,0.0078;-0.0092,0.999,-0.006;-0.007693,-0.00576,0.999];
        A_0=A_B_LVLH*A_LVLH_ECI_0;
        A_slew_B_LVLH=zeros(3);             % unused in this maneuver

    case 3
        omega_sc_0=[0,0,n];
        A_B_LVLH=eye(3);
        A_0=A_B_LVLH*A_LVLH_ECI_0;
        A_slew_B_LVLH=rot_matrix(deg2rad(-15.3),1)*rot_matrix(deg2rad(0),2)*rot_matrix(deg2rad(0),3)*eye(3);             % unused in this maneuver
end

% Reaction Wheels angular velocity in their frame, initial condition

om_rw_0 = [0,0,0];                           % [rad/s]

%% CONTROL

% State space rapresentation of the linearized system along the equilibrium
% configuration fixed to LVLH frame

A=[0,  (1-K_yaw)*n, 0, -K_yaw*n^2,  0,   0
   (K_roll-1)*n, 0, 0, 0, -4*K_roll*n^2, 0
   0,    0,      0, 0, 0,   -3*K_pitch*n^2    
   1,    0,      0, 0,                0, 0
   0,    1,      0, 0,                0, 0 
   0,    0,      1, 0,                0, 0];

B=[1/I(1,1),0,0
   0,1/I(2,2),0
   0,0,1/I(3,3)
       zeros(3)];

C=[0,0,0,1,0,0
   0,0,0,0,1,0
   0,0,0,0,0,1];

D=zeros(3,3);

% Optimal control weighting coefficients
% Linear quadratic regulator
% QQ=[[3000,0,0;0,10000,0;0,0,1000],zeros(3);zeros(3),[0.4,0,0;0,15,0;0,0,10]];
% RR=[15000,0,0;0,300000,0;0,0,1000];
% [K,S,P]=lqr(A,-B,QQ,RR);
% 
% % State observer 
% L=place(A',C',1.5*P)';


switch maneuver
    case 1
        % QQ=[[3000,0,0;0,10000,0;0,0,1000],zeros(3);zeros(3),[0.4,0,0;0,15,0;0,0,10]];
        % RR=[15000,0,0;0,300000,0;0,0,1000];
        % Poles of L: 1.5 poles of K

        K = [  -0.4671    0.0022    0.0000   -0.0052    0.0001   -0.0000
                0.0004   -0.2012    0.0000   -0.0000   -0.0071   -0.0000
               -0.0000   -0.0000   -1.1514   -0.0000   -0.0000   -0.1000];

        L = [  0.2157   -0.1255   -0.0001
              -0.1318    0.2078    0.0001
              -0.0000    0.0000    0.0010
               1.0261   -0.2800   -0.0001
              -0.2990    0.9541    0.0002
              -0.0001    0.0000    0.0758];

    case 2
        % QQ=[1*eye(3),zeros(3);zeros(3),1*eye(3)];
        % RR=100000*eye(3);
        % P_observer=3*P_gain
        % P.N.: Doesn't work with desaturation, one should switch to case 1

        K = [ -0.1081    0.0000   -0.0000   -0.0033    0.0000   -0.0000
               0.0000   -0.0580   -0.0000   -0.0000   -0.0033   -0.0000
              -0.0000    0.0000   -0.1041    0.0000    0.0000   -0.0033];

        L = [0.0317   -0.0001    0.0008
            -0.0007    0.0322    0.0009
             0.0008    0.0008    0.0181
             0.2641   -0.0791   -0.0047
             0.0759    0.2653    0.0072
             0.0132    0.0038    0.1908];

    case 3
        % QQ=[1500*eye(3),zeros(3);zeros(3),1.1*eye(3)];
        % RR=10000*eye(3);
        % P_observer=1.8*P_gain

        K = [ -0.4323    0.0000   -0.0000   -0.0105    0.0000   -0.0000
               0.0000   -0.4007    0.0000   -0.0000   -0.0105   -0.0000
              -0.0000   -0.0000   -0.4292    0.0000   -0.0000   -0.0105];

        L = [  0.0663    0.0001    0.0068
              -0.0019    0.0193   -0.0002
               0.0068   -0.0000    0.0219
               1.4076   -0.0003    0.1386
              -0.0024    0.4421   -0.0003
               0.1387   -0.0002    0.4949];
    end

if (executionType == 1)
    %% SIMULATION
    disp('Running simulink model...');
    % Params
    gamma_0 = [0;0;0];
    simfile = 'SAD_PROJECT.slx';
    %sim_time = 2*pi/n;
    %max_dt = .01; %[s] 'MaxStep', num2str(max_dt)
    % abs_tol = 1e-7;
    % rel_tol = 1e-7;
    % 
    % if bdIsLoaded(simfile)==1 %Checks if model is open so we can edit these params
    %     set_param(simfile, 'StopTime', 'sim_time',...
    %     'AbsTol', num2str(abs_tol), 'RelTol', num2str(rel_tol))
    % end
    
    % Execution
    out = sim(simfile);
    time = out.tout;
    gamma = out.gamma;

    disp('Postprocessing...');
    groundTrack = out.gamma;
    Rgt = 6380;
    for i=1:length(gamma)
        if i>length(gamma)/2
            gamma(i,1) = 2*pi-gamma(i,1);
        end
        groundTrack(i,:) = normalize(out.RR_Esc(i,1:3)','norm')'*Rgt;
    end
    disp('Plotting graphs...');

    % Plots
    simplePlot(time,{out.RR_Esc(:,1),out.RR_Esc(:,2),out.RR_Esc(:,3)},'Time [s]','r [km]','Earth-Spacecraft position',{'x','y','z'});
    simplePlot(time,{out.VV_Esc(:,1),out.VV_Esc(:,2),out.VV_Esc(:,3)},'Time [s]','r [km/s]','Earth-Spacecraft velocity',{'x','y','z'});
    simplePlot(time,{out.RR_ES(:,1),out.RR_ES(:,2),out.RR_ES(:,3)},'Time [s]','r [km]','Earth-Sun position',{'x','y','z'});
    simplePlot(time,{out.VV_ES(:,1),out.VV_ES(:,2),out.VV_ES(:,3)},'Time [s]','r [km/s]','Earth-Sun velocity',{'x','y','z'});
    simplePlot(time,{out.alpha_x,out.alpha_z},'Time [s]','Magnitude [deg]','Alpha vs time',{'x','z'});
    simplePlot(time,{out.alpha_x},'Time [s]','Yaw angle [deg]','Pointing error',{''});
    simplePlot(time,{out.ext_perturbations(:,1),out.ext_perturbations(:,2),out.ext_perturbations(:,3)},'Time [s]','F [N]','External perturbations',{'x','y','z'});
    simplePlot(time,{out.omega_rw(:,1),out.omega_rw(:,2),out.omega_rw(:,3)},'Time [s]','Magnitude [rad/s]','Reaction wheel Angular velocity',{'x','y','z'});
    simplePlot(time,{out.omega_sc(:,1),out.omega_sc(:,2),out.omega_sc(:,3)},'Time [s]','Magnitude [rad/s]','Spacecraft Angular velocity',{'wx','wy','wz'});
    simplePlot(time,{out.rw_perturbations.Data(:,1),out.rw_perturbations.Data(:,2),out.rw_perturbations.Data(:,3)},'Time [s]','H [Nm]','Reaction wheel perturbations',{'x','y','z'});
    simplePlot(time,{gamma(:,1),gamma(:,2),gamma(:,3)},'Time [s]','Magnitude []','Pointing direction',{'gammax','gammay','gammaz'});
    
end

%%
earthAngle = 0.00004667/180*pi*time;

%
step = 2;
framerate = 0.5;

tInterp = (min(time):step:max(time))';
rx = interp1(time, out.RR_Esc(:,1), tInterp, 'linear');
ry = interp1(time, out.RR_Esc(:,2), tInterp, 'linear');
rz = interp1(time, out.RR_Esc(:,3), tInterp, 'linear');


angx = interp1(time, gamma(:,1), tInterp, 'linear');
angy = interp1(time, gamma(:,2), tInterp, 'linear');
angz = interp1(time, gamma(:,3), tInterp, 'linear');

GTx = interp1(time, groundTrack(:,1), tInterp, 'linear');
GTy = interp1(time, groundTrack(:,2), tInterp, 'linear');
GTz = interp1(time, groundTrack(:,3), tInterp, 'linear');

earthAngle = interp1(time, earthAngle, tInterp, 'linear');
tInterp = round(tInterp*framerate);

data = table(tInterp,rx,ry,rz,angx,angy,angz,GTx,GTy,GTz,earthAngle,'VariableNames',{'t','x','y','z','angx','angy','angz','GTx','GTy','GTz','earthAngle'});
writetable(data,'trajectoryData.csv');

%% NESTED FUNCTIONS

function [rr,vv]=koe2rv(a,e,i,OM,om,theta,mu,p)

% FUNCTION NAME
%     koe2rv.m
% 
% DESCRIPTION
%     It computes the postion and velocity vectors given a set of orbital
%     elements. It accepts n multiple sets.
% 
% INPUT
%     MANDATORY:
%     in1 [1xn] - semimajor axis (km)
%     in2 [1xn] - eccentricity value (-)
%     in3 [1xn] - inclination (rad)
%     in4 [1xn] - right ascension of the ascending node (rad)
%     in5 [1xn] - argument of the perigee (rad)
%     in6 [1xn] - true anomaly (rad)
%     OPTIONAL:
%     in7 [1xn] - planetary constant (km^3/s^2), if not specified mu_Earth
%     in8 [1xm] - semilatus rectum for parabolic orbits (km)
%
% OUTPUT 
%     out1 - [3xn] position vector in GCFR (km)
%     out2 - [3xn] velocity vector in GCFR (km/s)
%
% INTERNAL PARAMETERS, Earth's planetary constant retrived via 
%     DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
%
% NESTED FUNCTION: rot_matrix
%
% CONTRIBUTORS
%     Andrea Binotto
% 
% VERSION
%     2023-10-22: First Version


% Initialization of p parameters for non parabolic orbits
    if nargin== 7
        p=[];

% Earth's planetary constant unless specified otherwise
    elseif nargin==6
        mu=3.98600433e+5;                                      %(km^3/s^2)
        p=[];
    end

% radius in perifocal frame
    rr_pf=a.*(1-e.^2)./(1+e.*cos(theta)).*[1;0;0];
   
    % parabolic orbit
    if max(isnan(a))==1
       rr_pf(:,isnan(a))=p./(1+e(isnan(a)).*cos(theta(isnan(a)))).*[1;0;0];    
    end

% velocity norm
    h=sqrt(mu*(a.*(1-e.^2)));
    h(isnan(a))=sqrt(mu*p);
    vv_pf=mu./h.*[e.*sin(theta);1+e.*cos(theta);zeros(1,length(a))];

% GCRF to PF rotation matrix
    R_OM=zeros(3,3,length(a));
    R_i=zeros(3,3,length(a));
    R_om_theta=zeros(3,3,length(a));
    rr=zeros(3,length(a));
    vv=zeros(3,length(a));

for j=1:length(a)
    R_OM(:,:,j)= rot_matrix(-OM(j),3);
    R_i(:,:,j)= rot_matrix(-i(j),1);
    R_om_theta(:,:,j)= rot_matrix(-om(j)-theta(j),3);
    
    % PF to GCRF
    rr(:,j)= R_OM(:,:,j)'*R_i(:,:,j)'*R_om_theta(:,:,j)'*rr_pf(:,j);
    vv(:,j)= R_OM(:,:,j)'*R_i(:,:,j)'*R_om_theta(:,:,j)'*vv_pf(:,j);
end


% NESTED FUNCTION: ROTATION MATRIX

    function [R]=rot_matrix(angle,direction)
    
    % FUNCTION NAME
    %     rot_matrix.m
    % 
    % DESCRIPTION
    %     It computes the matrix for a counter clockwise rotation of amplitude 
    %     "angle" along principal "direction".
    % 
    % INPUT
    %     in1 [1] - Rotation angle, counter clockwise (rad).
    %     in2 [1] - Axis of rotation, either (1=x,2=y,3=z) .
    % 
    % OUTPUT 
    %     out1 [3x3] - Rotation matrix.
    %
    % CONTRIBUTORS
    %     Andrea Binotto
    % 
    % VERSION
    %     2023-10-22: First Version
    
    switch direction
        case 1
            R=[1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)];
        case 2
            R=[cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
        case 3
            R=[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
    end
    
    end


end

function [dydt] = twobody_J2 (t,s,mu)

% FUNCTION NAME
%     2bp_prt
% 
% DESCRIPTION
%     It computes restricted two body problem with chosen perturbations
% 
% INPUT
%     in1 [1] - Time (can be omitted as the system is autonomous)
%     in2 [6x1] - State vector [position vector ; velocity vector]
%     in3 [1] - Planet's gravitational constant (km^3/s^2)
% 
% OUTPUT 
%     out1 - [6x1] Differential problem of [r;v]
%
% INTERNAL PARAMETERS, Earth's planetary constant retrived via 
%     DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
%
% CONTRIBUTORS
%     Andrea Binotto
% 
% VERSION
%     2023-10-15: First Version

% Earth's planetary constant unless specified otherwise
    if nargin== 2 
       mu=3.98600433e+5;                                    %(km^3/s^2)
    end

% J2 perturbation

J2=1.082626925638815e-3;     %Gravitatonal field constant of the Earth
R=6378.1363;                 %Earth's equatorial radius [km]

% a_J2=1.5*J2*mu*R^2/(norm(s(1:3))^4)*[
% s(1)/norm(s(1:3))*(5*s(3)^2/(norm(s(1:3))^2)-1);
% s(2)/norm(s(1:3))*(5*s(3)^2/(norm(s(1:3))^2)-1);
% s(3)/norm(s(1:3))*(5*s(3)^2/(norm(s(1:3))^2)-3)];


a_J2=1.5*J2*mu*R^2/(norm(s(1:3))^7)*[
s(1)*(5*s(3)^2-(norm(s(1:3))^2));
s(2)*(5*s(3)^2-(norm(s(1:3))^2));
s(3)*(5*s(3)^2-3*(norm(s(1:3))^2))];

% Keplerian Motion
dr=s(4:6);
dv=-mu/(norm(s(1:3))^3)*s(1:3);

% Pertubrated equations of motion
dydt=[dr;dv+a_J2];

end

function [theta_e] = orbit_propagation(a,e,theta_s,delta_t,mu,p)

% FUNCTION NAME
%     orbit_propagation.m
%
% DESCRIPTION
%     It computes the true anomaly reached along an orbit given the time
%     of flight and the initial position. It solves Kepler problem by a
%     Newton-Raphson method.
%
% INPUT
%     MANDATORY:
%     in1 [1] - semimajor axis (km)
%     in2 [1] - eccentricity value (-)
%     in3 [1] - true anomaly of starting point (rad)
%     in4 [1] - time of flight (s)
%     OPTIONAL:
%     in5 [1] - planetary constant (km^3/s^2), if not specified mu_Earth
%     in6 [1] - semilatus rectum for parabolic orbits only
%
% OUTPUT
%     out1 [1] - true anomaly after given time.
%
% NESTED FUNCTION
%     newton_method.m
%
% CONTRIBUTORS
%     Andrea Binotto
%
% VERSION
%     2023-10-26: First Version
% Different for kind of orbit



% Initialization of p parameters for non parabolic orbits
    if nargin== 5
        p=[];

% Earth's planetary constant unless specified otherwise
    elseif nargin==4
        mu=3.98600433e+5;                           % (km^3/s^2)
        p=[];
    end


% orbit type discrimination

if e<1                                              % ecliptic orbit

    E_s=2*atan(sqrt((1-e)/(1+e))*tan(theta_s/2));   % eccentric anomaly
    n_e=sqrt(mu/(a^3));                             % mean motion (rad/s)
    T=2*pi/n_e;

    rev=0;                                          % n. of full orbit
    while delta_t>T
        delta_t=delta_t-T;
        rev=rev+1;
    end

    M_e=n_e*delta_t;                                % mean anomaly (rad)
    f=@(E) E-e*sin(E)-M_e-E_s+e*sin(E_s);           % function
    df=@(E) 1-e*cos(E);                             % function derivative
    toll=1e-15;                            % tolerance for absolute error

    % choice of initial guess, as suggested in (Prussing and Conway, 2013)
    if mod(M_e,2*pi)>pi             % from apocenter to pericenter branch
        E_e0=M_e-e/2;                 % E behind M
    else                            % from pericenter to apocenter branch
        E_e0=M_e+e/2;                 % E beyond M
    end

    E_e=newton_method(E_e0,f,df,toll);           % eccentric anomaly (rad)
    theta_e=mod(2*atan(sqrt((1+e)/(1-e))*tan(E_e/2)),2*pi); % true anomaly
    
    % progressive orbit motion
    if theta_e<theta_s
        rev=rev+1;
    end
    
    theta_e=theta_e+rev*2*pi;


elseif e==1        % parabolic orbit
    
    n_p=2*sqrt(mu/(p^3));                           % mean motion (rad/s)
    M_p=n_p*delta_t;                                % mean anomaly (rad)
    D_s=tan(theta_s/2);                             % parabolic anomaly
    f=@(D) 0.5*(D+D^3/3)-0.5*(D_s+D_s^3/3)-M_p;     % function
    D_e=fzero(f,0);                                 % parabolic anomaly
    
    theta_e=2*atan(D_e);


else               % hyperbolic orbit

    % belonging to physical branch check
    theta_inf=acos(-1/e);
    if abs(theta_s)>theta_inf
        error 'point belongs to not-physical branch'
    end

    n_h=sqrt(mu/(abs(a)^3));                       % mean motion (rad/s)
    M_h=n_h*delta_t;                               % mean anomaly (rad)
    F_s=2*atanh(sqrt((e-1)/(e+1))*tan(theta_s/2)); % hyperbolic anomaly
    f=@(F) -F+e*sinh(F)+F_s-e*sinh(F_s)-M_h;
    df=@(F) -1+e*cosh(F);
    toll=1e-13;                            % tolerance for absolute error

    % choice of initial guess, as suggested in (Prussing and Conway, 2013)
    if mod(M_h,2*pi)>pi             
        F_e0=M_h-e/2;               
    else                            
        F_e0=M_h+e/2;                 
    end

    F_e=newton_method(F_e0,f,df,toll);          % hyperbolic anomaly (rad)
    theta_e=2*atan(sqrt((1+e)/(e-1))*tanh(F_e/2)); % true anomaly
end

function [x]=newton_method(x0,fun,dfun,toll, mol)

% FUNCTION NAME
%     newton_method.m
%
% DESCRIPTION
%     Newton's method for root finding of a function. Stopping criterion
%     based on absolute error.
%
% INPUT
%     MANDATORY:
%     in1 [1] - initial guess
%     in2 [1] - function (anonymus)
%     in3 [1] - function first derivative (anonymus)
%     in4 [1] - tolerance
%     OPTIONAL:
%     in5 [1] - moltiplicity of root

% OUTPUT
%     out1 [1] - root.
%
% CONTRIBUTORS
%     Andrea Binotto
%
% VERSION
%     2023-10-26: First Version


% single moltiplicity unless specified otherwise 

if nargin == 4
    mol=1;
end

abs_err=toll+1;                                     % absolute tolerance
it=0;                                               % n. of iteration
x_old= x0;

while abs_err>toll                                  % stopping criterion
   df_eval=dfun(x_old);
   if df_eval==0
      error(' flat function');
   else
      x=x_old-mol*fun(x_old)./df_eval;               % newton's algorithm
      abs_err=abs(x-x_old);                          
      x_old= x;
      it=it+1;
   end
end

if it>100
    disp ('it converges after more than 100 iterations')
end
end
end

function [R]=rot_matrix(angle,direction)

% FUNCTION NAME
%     rot_matrix.m
% 
% DESCRIPTION
%     It computes the matrix for a counter clockwise rotation of amplitude 
%     "angle" along principal "direction".
% 
% INPUT
%     in1 [1] - Rotation angle, counter clockwise (rad).
%     in2 [1] - Axis of rotation, either (1=x,2=y,3=z) .
% 
% OUTPUT 
%     out1 [3x3] - Rotation matrix.
%
% CONTRIBUTORS
%     Andrea Binotto
% 
% VERSION
%     2023-10-22: First Version

switch direction
    case 1
        R=[1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)];
    case 2
        R=[cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
    case 3
        R=[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
end
end