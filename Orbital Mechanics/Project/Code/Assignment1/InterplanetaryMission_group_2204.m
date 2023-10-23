% Main script for Orbital Mechanics project
% Academic Year: 2022/2023
% Group ID: 2204

clc
close all
clear
%% Data
G = astroConstants(1);
Sun = struct;
Sun.mu = astroConstants(4);
Sun.m = Sun.mu / G;

% For ode solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Earliest departure from Earth and latest arrival to NEO 85
earliest_dep = [2026, 9, 30, 0, 0, 0];
latest_arr = [2061, 3, 28, 23, 59, 59];
mjd2000_earliest_dep = date2mjd2000(earliest_dep); % [d]
mjd2000_latest_arr = date2mjd2000(latest_arr); % [d]

% Departure planet: Earth
Earth = struct;
Earth.ID = 3;
[kep_E, ~] = uplanet(mjd2000_earliest_dep, Earth.ID);
Earth.T = 2*pi * sqrt(kep_E(1)^3/Sun.mu);
Earth.n = 2*pi / Earth.T;
clear kep_E
% Fly-by planet: Saturn
Saturn = struct;
Saturn.ID = 6;
Saturn.mu = astroConstants(16);
Saturn.m = Saturn.mu / G;
Saturn.R = astroConstants(26);
Saturn.h_atm = 59.5; % Saturn atmosphere scale height
Saturn.h_outring = 480000 - Saturn.R;
[kep_S, ~] = uplanet(mjd2000_earliest_dep, Saturn.ID);
Saturn.T = 2*pi * sqrt(kep_S(1)^3/Sun.mu);
Saturn.n = 2*pi / Saturn.T;
Saturn.r_SOI = kep_S(1) * (Saturn.m/Sun.m)^(2/5);
clear kep_S
% Arrival NEO: 85
NEO = struct;
NEO.ID = 85;
[kep_NEO, ~, ~] = ephNEO(mjd2000_earliest_dep,85);
NEO.T = 2*pi * sqrt(kep_NEO(1)^3/Sun.mu);
NEO.n = 2*pi / NEO.T;
clear kep_NEO

% Other data
excess_v = 100;

%% 3-body configuration repetition
% Sinodic years
Tsyn_ES = 2*pi / abs(Earth.n-Saturn.n); % [s]
Tsyn_SNEO = 2*pi / abs(NEO.n-Saturn.n); % [s]

rep = 1:20;
Tsyn_ES_vec =  seconds2days(Tsyn_ES) .* rep; % Tsyn between E and S propagated in time until its 19th repetition
Tsyn_SNEO_vec = seconds2days(Tsyn_SNEO) .* rep; % Tsyn between S and NEO propagated in time until its 19th repetition

diff_act = 30; i_act = 0; j_act =0; k = 0; % Initilising variables
for i = 1:length(rep)
    for j = 1:length(rep)
        diff = abs(Tsyn_SNEO_vec(j) - Tsyn_ES_vec(i));
        if diff < 30 && diff < diff_act % Difference must be less than 30 days
            diff_act = diff;
            i_act = i;
            j_act = j;
        end
    end
end
if diff_act == 30
    error('Total synodic period not found.')
else
    Tsyn_tot = Tsyn_ES * rep(i_act); % [s]
end

%% Time windows
window_dep_mjd = [mjd2000_earliest_dep, mjd2000_earliest_dep + seconds2days(Tsyn_tot)];
window_fb_mjd = [window_dep_mjd(1)+1, mjd2000_latest_arr-1];
window_arr_mjd = [window_fb_mjd(1)+1, mjd2000_latest_arr];
dep_window = linspace(window_dep_mjd(1), window_dep_mjd(2), 50);
fb_window = linspace(window_fb_mjd(1), window_fb_mjd(2), 50);
arr_window = linspace(window_arr_mjd(1), window_arr_mjd(2), 50);
% Search for minimum parabolic flight time
tpar_min = min_parabolic_time(dep_window,fb_window,arr_window,Earth.ID,Saturn.ID,NEO.ID,Sun.mu);
% Definition of time windows
window_dep_mjd = [mjd2000_earliest_dep, mjd2000_earliest_dep + seconds2days(Tsyn_tot)];
window_fb_mjd = [window_dep_mjd(1) + seconds2days(tpar_min), mjd2000_latest_arr - seconds2days(tpar_min)];
window_arr_mjd = [window_fb_mjd(1) + seconds2days(tpar_min), mjd2000_latest_arr];

dep_window = linspace(window_dep_mjd(1), window_dep_mjd(2), 81);
fb_window = linspace(window_fb_mjd(1), window_fb_mjd(2), 35);
arr_window = linspace(window_arr_mjd(1), window_arr_mjd(2), 125);

%% Mission (No constraints on maximum departure date)
IDs = [Earth.ID; Saturn.ID; NEO.ID];
opt_lambert_solv = struct;
opt_lambert_solv.orbitType1 = 0; opt_lambert_solv.orbitType2 = 0;
opt_lambert_solv.RegionsNumber = 3; opt_lambert_solv.tpar_min = tpar_min;
% Computing minimum cost and dates for departure, fly-by and arrival
[dv_min, date_mjd2000] = transfer_design_3B(dep_window,fb_window,arr_window,...
    IDs,excess_v,Saturn.h_outring,opt_lambert_solv);

%% Data computation
% Dates
date_sec = zeros(3,1); date = zeros(3,6);
for i = 1:3
    date_sec(i) = days2seconds(date_mjd2000(i));
    date(i,:) = mjd20002date(date_mjd2000(i));
end
% Keplerian and cartesian parameters of space bodies at departure
[kep_E, ~] = uplanet(date_mjd2000(1), Earth.ID); 
[kep_S, ~] = uplanet(date_mjd2000(2), Saturn.ID);
[kep_NEO, ~, ~] = ephNEO(date_mjd2000(3), NEO.ID);
[r_E, v_E] = kep2car(kep_E, Sun.mu); 
[r_S, v_S] = kep2car(kep_S, Sun.mu); 
[r_NEO, v_NEO] = kep2car(kep_NEO, Sun.mu);
% Time of flights
ToF1 = date_sec(2) - date_sec(1); ToF2 = date_sec(3) - date_sec(2);
ToF_tot = ToF1 + ToF2;
% Lambert solver
[~, ~, ~, ~, v1_t1, v2_t1, tpar1, ~] = lambertMR(r_E, r_S, ToF1, Sun.mu, opt_lambert_solv.orbitType1, ...
                                        0, 0, 2);
[~, ~, ~, ~, v1_t2, v2_t2, tpar2, ~] = lambertMR(r_S, r_NEO, ToF2, Sun.mu, opt_lambert_solv.orbitType2, ...
                                        0, 0, 2);
v1_t1 = v1_t1'; v2_t1 = v2_t1'; v1_t2 = v1_t2'; v2_t2 = v2_t2';
kep_t1_i = car2kep(r_E, v1_t1, Sun.mu); kep_t1_f = car2kep(r_S, v2_t1, Sun.mu);
kep_t2_i = car2kep(r_S, v1_t2, Sun.mu); kep_t2_f = car2kep(r_NEO, v2_t2, Sun.mu);
% Costs and flyby
deltaV1 = abs(norm(v1_t1 - v_E)); deltaV2 = abs(norm(v_NEO - v2_t2));
deltaV_fb = v1_t2 - v2_t1; % Change of heliocentric velocity thanks to flyby
v_inf_m = v2_t1 - v_S; v_inf_p = v1_t2 - v_S; % s/c velocities relative to the planet
[r_p, delta, deltaV_ga, kep_hyp1, kep_hyp2, centres] = power_gravity_assist(v_inf_m, v_inf_p, Saturn.ID, Saturn.h_outring);
% Hyperbolae parameters
u1 = cross(r_p, v_inf_m) / norm(cross(r_p, v_inf_m));
v_p1_unitvec = cross(u1, r_p) / norm(cross(u1, r_p)); % Unitvector of velocity @pericentre hyp1
v_p2_unitvec = cross(u1, r_p) / norm(cross(u1, r_p)); % Unitvector of velocity @pericentre hyp1
h_hyp1 = sqrt(Saturn.mu * kep_hyp1(1) * (1-kep_hyp1(2)^2));
h_hyp2 = sqrt(Saturn.mu * kep_hyp2(1) * (1-kep_hyp2(2)^2));
v_p1 = Saturn.mu / h_hyp1 * (1 + kep_hyp1(2)) * v_p1_unitvec; % Velocity @pericentre hyp1
v_p2 = Saturn.mu / h_hyp2 * (1 + kep_hyp2(2)) * v_p2_unitvec; % Velocity @pericentre hyp2
kep_hyp1 = car2kep(r_p, v_p1, Saturn.mu);
kep_hyp2 = car2kep(r_p, v_p2, Saturn.mu);
% Time of Flight flyby
th_hyp = @(a,e) acos(1/e * (a*(1-e^2)/Saturn.r_SOI - 1)); % True anomaly
E1_bar = @(a,e) 2*atanh(tan(th_hyp(a,e)/2) * sqrt((e-1)/(1+e))); % Eccentric anomaly
deltaT = @(a,e) (e*sinh(E1_bar(a,e)) - E1_bar(a,e)) * sqrt((-a)^3/Saturn.mu); % ToF
deltaT1_hyp = deltaT(kep_hyp1(1), kep_hyp1(2)); deltaT2_hyp = deltaT(kep_hyp2(1), kep_hyp2(2));
ToF_flyby = deltaT1_hyp + deltaT2_hyp; ToF_flyby_hours = ToF_flyby / (60*60);

%% Mission (with constraint on maximum departure date)
window_dep_mjd = [mjd2000_earliest_dep, mjd2000_earliest_dep + (365.25*10)];
window_fb_mjd = [window_dep_mjd(1) + seconds2days(tpar_min), mjd2000_latest_arr - seconds2days(tpar_min)];
window_arr_mjd = [window_fb_mjd(1) + seconds2days(tpar_min), mjd2000_latest_arr];

dep_window = linspace(window_dep_mjd(1), window_dep_mjd(2), 61);
fb_window = linspace(window_fb_mjd(1), window_fb_mjd(2), 35);
arr_window = linspace(window_arr_mjd(1), window_arr_mjd(2), 125);
opt_lambert_solv.RegionsNumber = 2;
% Computing minimum cost and dates for departure, fly-by and arrival
[dv_min_opt2, date_mjd2000_opt2] = transfer_design_3B(dep_window,fb_window,arr_window,IDs,excess_v,Saturn.h_outring,opt_lambert_solv);

%% Data computation
% Dates
date_sec_opt2 = zeros(3,1); date_opt2 = zeros(3,6);
for i = 1:3
    date_sec_opt2(i) = days2seconds(date_mjd2000_opt2(i));
    date_opt2(i,:) = mjd20002date(date_mjd2000_opt2(i));
end
% Keplerian and cartesian parameters of space bodies at departure
[kep_E_opt2, ~] = uplanet(date_mjd2000_opt2(1), Earth.ID); 
[kep_S_opt2, ~] = uplanet(date_mjd2000_opt2(2), Saturn.ID);
[kep_NEO_opt2, ~, ~] = ephNEO(date_mjd2000_opt2(3), NEO.ID);
[r_E_opt2, v_E_opt2] = kep2car(kep_E_opt2, Sun.mu); 
[r_S_opt2, v_S_opt2] = kep2car(kep_S_opt2, Sun.mu); 
[r_NEO_opt2, v_NEO_opt2] = kep2car(kep_NEO_opt2, Sun.mu);
% Time of flights
ToF1_opt2 = date_sec_opt2(2) - date_sec_opt2(1); ToF2_opt2 = date_sec_opt2(3) - date_sec_opt2(2);
ToF_tot_opt2 = ToF1_opt2 + ToF2_opt2;
% Lambert solver
[~, ~, ~, ~, v1_t1_opt2, v2_t1_opt2, tpar1_opt2, ~] = lambertMR(r_E_opt2, r_S_opt2, ToF1_opt2, Sun.mu, ...
                                opt_lambert_solv.orbitType1, 0, 0, 2);
[~, ~, ~, ~, v1_t2_opt2, v2_t2_opt2, tpar2_opt2, ~] = lambertMR(r_S, r_NEO, ToF2, Sun.mu, ...
                                opt_lambert_solv.orbitType2, 0, 0, 2);
v1_t1_opt2 = v1_t1_opt2'; v2_t1_opt2 = v2_t1_opt2'; 
v1_t2_opt2 = v1_t2_opt2'; v2_t2_opt2 = v2_t2_opt2';
kep_t1_opt2_i = car2kep(r_E_opt2, v1_t1_opt2, Sun.mu); kep_t1_opt2_f = car2kep(r_S_opt2, v2_t1_opt2, Sun.mu);
kep_t2_opt2_i = car2kep(r_S_opt2, v1_t2_opt2, Sun.mu); kep_t2_opt2_f = car2kep(r_NEO_opt2, v2_t2_opt2, Sun.mu);
% Costs and flyby
deltaV1_opt2 = abs(norm(v1_t1_opt2 - v_E_opt2)); deltaV2_opt2 = abs(norm(v_NEO_opt2 - v2_t2_opt2));
deltaV_fb_opt2 = v1_t2_opt2 - v2_t1_opt2; % Change of heliocentric velocity thanks to flyby
v_inf_m_opt2 = v2_t1_opt2 - v_S_opt2; v_inf_p_opt2 = v1_t2_opt2 - v_S_opt2; % s/c velocities relative to the planet
[r_p_opt2, delta_opt2, deltaV_ga_opt2, ~, ~, centres_opt2] = ...
    power_gravity_assist(v_inf_m_opt2, v_inf_p_opt2, Saturn.ID, Saturn.h_outring);
%% Velocities triangle
vS_vec = [[0;0;0], v_S];
v2_t1_vec = [[0;0;0], v2_t1];
v1_t2_vec = [[0;0;0], v1_t2];
v_inf_m_vec = [[0;0;0], v_inf_m] + v_S;
v_inf_p_vec = [[0;0;0], v_inf_p] + v_S;
figure()
hold on
grid on
axis equal
plot3(vS_vec(1,:),vS_vec(2,:),vS_vec(3,:),'Color', "#EDB120",'LineWidth',2)
plot3(v2_t1_vec(1,:),v2_t1_vec(2,:),v2_t1_vec(3,:),'Color', "#77AC30",'LineWidth',2)
plot3(v_inf_m_vec(1,:),v_inf_m_vec(2,:),v_inf_m_vec(3,:),'r')
plot3(v1_t2_vec(1,:),v1_t2_vec(2,:),v1_t2_vec(3,:),'Color', "#4DBEEE",'LineWidth',2)
plot3(v_inf_p_vec(1,:),v_inf_p_vec(2,:),v_inf_p_vec(3,:),'g')
legend('$\underline{V}_S$','$\underline{V}^-$','$\underline{v}^-_{\infty}$',...
    '$\underline{V}^+$','$\underline{v}^+_{\infty}$','Interpreter','latex')
xlabel('$x$ component [km/s]','Interpreter','latex');
ylabel('$y$ component [km/s]','Interpreter','latex');
zlabel('$z$ component [km/s]','Interpreter','latex');
title('Velocities triangle')

%% Data for plots
tplot_E = linspace(date_mjd2000(1), date_mjd2000(1) + seconds2days(Earth.T), 5000);
tplot_S = linspace(date_mjd2000(1), date_mjd2000(1) + seconds2days(Saturn.T), 50000);
tplot_NEO = linspace(date_mjd2000(1), date_mjd2000(1) + seconds2days(NEO.T), 5000);
tplot_t1 = linspace(date_sec(1), date_sec(2), 10000);
tplot_t2 = linspace(date_sec(2), date_sec(3), 10000);
y0_t1 = [r_E; v1_t1]; y0_t2 = [r_S; v1_t2];

Y_E = zeros(length(tplot_E), 3); Y_S = zeros(length(tplot_S), 3); Y_NEO = zeros(length(tplot_NEO), 3);
for i = 1:length(tplot_E)
    [kep_E_i, ~] = uplanet(tplot_E(i), Earth.ID);
    [Y_E(i,:), ~] = kep2car(kep_E_i, Sun.mu);
end
for i = 1:length(tplot_S)
    [kep_S_i, ~] = uplanet(tplot_S(i), Saturn.ID);
    [Y_S(i,:), ~] = kep2car(kep_S_i, Sun.mu);
end
for i = 1:length(tplot_NEO)
    [kep_NEO_i, ~, ~] = ephNEO(tplot_NEO(i), NEO.ID);
    [Y_NEO(i,:), ~] = kep2car(kep_NEO_i, Sun.mu);
end

[~, Y_t1] = ode113(@(t,y) ode_2bp(t,y,Sun.mu), tplot_t1, y0_t1, options);
[~, Y_t2] = ode113(@(t,y) ode_2bp(t,y,Sun.mu), tplot_t2, y0_t2, options);

%% Plot heliocentric trajectories
% Derparture
figure()
hold on
grid on
plot3(Y_E(:,1), Y_E(:,2), Y_E(:,3),'Color', "#0072BD", "LineWidth", 1.5)
plot3(Y_S(:,1), Y_S(:,2), Y_S(:,3),'Color', "#EDB120", "LineWidth", 1.5)
plot3(Y_NEO(:,1), Y_NEO(:,2), Y_NEO(:,3),'Color', "#7E2F8E", "LineWidth", 1.5)
plot3(Y_t1(:,1), Y_t1(:,2), Y_t1(:,3),'--','Color', "#77AC30")
plot3(Y_t2(:,1), Y_t2(:,2), Y_t2(:,3),'--','Color', "#4DBEEE")
celestial_body(11,0,0,0,75)
celestial_body(3,Y_E(1,1),Y_E(1,2),Y_E(1,3), 5000)
celestial_body(6,Y_S(1,1),Y_S(1,2),Y_S(1,3), 800)
scatter3(Y_NEO(1,1),Y_NEO(1,2),Y_NEO(1,3),[],[0.4940 0.1840 0.5560],"filled")
legend('Earth orbit', 'Saturn orbit', 'NEO orbit','First transfer trajectory', ...
    'Second transfer trajectory','','','','','Location','best')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
% title('Departure')
% Fly-by
[kepE_fb,~] = uplanet(date_mjd2000(2),Earth.ID); [rE_fb, ~] = kep2car(kepE_fb,Sun.mu);
[kepS_fb,~] = uplanet(date_mjd2000(2),Saturn.ID); [rS_fb, ~] = kep2car(kepS_fb,Sun.mu);
[kepNEO_fb,~,~] = ephNEO(date_mjd2000(2),NEO.ID); [rNEO_fb, ~] = kep2car(kepNEO_fb,Sun.mu);
figure()
hold on
grid on
plot3(Y_E(:,1), Y_E(:,2), Y_E(:,3),'Color', "#0072BD", "LineWidth", 1.5)
plot3(Y_S(:,1), Y_S(:,2), Y_S(:,3),'Color', "#EDB120", "LineWidth", 1.5)
plot3(Y_NEO(:,1), Y_NEO(:,2), Y_NEO(:,3),'Color', "#7E2F8E", "LineWidth", 1.5)
plot3(Y_t1(:,1), Y_t1(:,2), Y_t1(:,3),'Color', "#77AC30", "LineWidth", 1.5)
plot3(Y_t2(:,1), Y_t2(:,2), Y_t2(:,3),'--','Color', "#4DBEEE")
celestial_body(11,0,0,0,75)
celestial_body(3,rE_fb(1),rE_fb(2),rE_fb(3), 5000)
celestial_body(6,rS_fb(1),rS_fb(2),rS_fb(3), 750)
scatter3(rNEO_fb(1),rNEO_fb(2),rNEO_fb(3),[],[0.4940 0.1840 0.5560],"filled")
legend('Earth orbit', 'Saturn orbit', 'NEO orbit','First transfer trajectory',...
    'Second transfer trajectory','','','','','Location','best')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
% title('Flyby')
% Arrival
[kepE_arr,~] = uplanet(date_mjd2000(3),Earth.ID); [rE_arr, ~] = kep2car(kepE_arr,Sun.mu);
[kepS_arr,~] = uplanet(date_mjd2000(3),Saturn.ID); [rS_arr, ~] = kep2car(kepS_arr,Sun.mu);
[kepNEO_arr,~,~] = ephNEO(date_mjd2000(3),NEO.ID); [rNEO_arr, ~] = kep2car(kepNEO_arr,Sun.mu);
figure()
hold on
grid on
plot3(Y_E(:,1), Y_E(:,2), Y_E(:,3),'Color', "#0072BD", "LineWidth", 1.5)
plot3(Y_S(:,1), Y_S(:,2), Y_S(:,3),'Color', "#EDB120", "LineWidth", 1.5)
plot3(Y_NEO(:,1), Y_NEO(:,2), Y_NEO(:,3),'Color', "#7E2F8E", "LineWidth", 1.5)
plot3(Y_t1(:,1), Y_t1(:,2), Y_t1(:,3),'Color', "#77AC30", "LineWidth", 1.5)
plot3(Y_t2(:,1), Y_t2(:,2), Y_t2(:,3),'Color', "#4DBEEE", "LineWidth", 1.5)
celestial_body(11,0,0,0,75)
celestial_body(3,rE_arr(1),rE_arr(2),rE_arr(3), 5000)
celestial_body(6,rS_arr(1),rS_arr(2),rS_arr(3), 750)
scatter3(rNEO_arr(1),rNEO_arr(2),rNEO_arr(3),[],[0.4940 0.1840 0.5560],"filled")
legend('Earth orbit', 'Saturn orbit', 'NEO orbit','First transfer trajectory',...
    'Second transfer trajectory','','','','','Location','best')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
% title('Arrival')

%% Hyperbola plots

% From Inertial to perifocal
i1 = kep_hyp1(3); OM1 = kep_hyp1(4); om1 = kep_hyp1(5);
i2 = kep_hyp2(3); OM2 = kep_hyp2(4); om2 = kep_hyp2(5);
ROM = [cos(OM1), sin(OM1), 0; -sin(OM1), cos(OM1), 0; 0, 0, 1];
Ri = [1, 0, 0; 0, cos(i1), sin(i1); 0, -sin(i1), cos(i1)];
Rom = [cos(om1), sin(om1), 0; -sin(om1), cos(om1), 0; 0, 0, 1];
R1 = Rom * Ri * ROM;
ROM = [cos(OM2), sin(OM2), 0; -sin(OM2), cos(OM2), 0; 0, 0, 1];
Ri = [1, 0, 0; 0, cos(i2), sin(i2); 0, -sin(i2), cos(i2)];
Rom = [cos(om2), sin(om2), 0; -sin(om2), cos(om2), 0; 0, 0, 1];
R2 = Rom * Ri * ROM;

% Asypmtotes
v_inf_m_unitvec = v_inf_m / norm(v_inf_m);
v_inf_p_unitvec = v_inf_p / norm(v_inf_p);
asyp1 = [centres(:,1) + 3*1e6*(-v_inf_m_unitvec), centres(:,1)]';
asyp2 = [centres(:,2) + 3*1e6*v_inf_p_unitvec, centres(:,2)]';
aps_line = [centres(:,2) + 5e1*(-r_p), centres(:,2)]';

y0_hyp1 = [r_p; -v_p1]; y0_hyp2 = [r_p; v_p2];
tspan_hyp = linspace(0, days2seconds(50), 1000);
[~, Y_hyp1] = ode113(@(t,y) ode_2bp(t,y,Saturn.mu), tspan_hyp, y0_hyp1, options);
[~, Y_hyp2] = ode113(@(t,y) ode_2bp(t,y,Saturn.mu), tspan_hyp, y0_hyp2, options);

% Inertial
figure()
grid on
opts_planet = struct;
opts_planet.RefPlane = 'ecliptic'; opts_planet.Units = 'km';
planet3D('Saturn',opts_planet);
hold on
p1 = plot3(Y_hyp1(:,1), Y_hyp1(:,2), Y_hyp1(:,3),'Color', "#77AC30", "LineWidth", 1.5);
p2 = plot3(Y_hyp2(:,1), Y_hyp2(:,2), Y_hyp2(:,3),'Color', "#4DBEEE", "LineWidth", 1.5);
p3 = plot3(asyp1(:,1), asyp1(:,2), asyp1(:,3),'--','Color', "#77AC30");
p4 = plot3(asyp2(:,1), asyp2(:,2), asyp2(:,3),'--','Color', "#4DBEEE");
p5 = plot3(aps_line(:,1), aps_line(:,2), aps_line(:,3),'k--');
plegend = [p1, p2, p3, p4, p5];
legend(plegend,'Incoming hyperbola','Outcoming hyperbola','Asymptote incoming hyperbola', ...
    'Asymptote outcoming hyperbola', 'Apse line','Location','best')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

for i = 1:length(tspan_hyp)
    Y_hyp1(i,1:3) = R1 * Y_hyp1(i,1:3)';
    Y_hyp2(i,1:3) = R2 * Y_hyp2(i,1:3)';
end
for i = 1:2
    asyp1(i,:) = R1 *  asyp1(i,:)';
    asyp2(i,:) = R2 *  asyp2(i,:)';
    aps_line(i,:) = R1 *  aps_line(i,:)';
end
% Perifocal
figure()
hold on
axis equal
% celestial_body(6)
opts_planet.RefPlane = 'equatorial';
planet_surface = planet3D('Saturn',opts_planet);
hold on
p1 = plot3(Y_hyp1(:,1), Y_hyp1(:,2), Y_hyp1(:,3),'Color', "#77AC30", "LineWidth", 1.5);
p2 = plot3(Y_hyp2(:,1), Y_hyp2(:,2), Y_hyp2(:,3),'Color', "#4DBEEE", "LineWidth", 1.5);
p3 = plot3(asyp1(:,1), asyp1(:,2), asyp1(:,3),'--','Color', "#77AC30");
p4 = plot3(asyp2(:,1), asyp2(:,2), asyp2(:,3),'--','Color', "#4DBEEE");
p5 = plot3(aps_line(:,1), aps_line(:,2), aps_line(:,3),'k--');
plegend = [p1, p2, p3, p4, p5];
grid on
legend(plegend,'Incoming hyperbola','Outcoming hyperbola','Asymptote incoming hyperbola', ...
    'Asymptote outcoming hyperbola', 'Apse line','Location','best')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

%% Porkchop plots
dep_window_ES = [mjd2000_earliest_dep; mjd2000_earliest_dep + seconds2days(Saturn.T)];
arr_window_ES = [mjd2000_earliest_dep; mjd2000_earliest_dep + seconds2days(Saturn.T)];
[dep_ES,arr_ES,dv_ES,delta_t_ES,r1_ES,v1_ES,~,~,dv1_ES,dv2_ES] = ...
    transfer_design(dep_window_ES,arr_window_ES,Earth.ID,Saturn.ID,Sun.mu,1);

dep_window_SNEO = [mjd2000_earliest_dep; mjd2000_earliest_dep + seconds2days(Saturn.T)];
arr_window_SNEO = [mjd2000_earliest_dep; mjd2000_earliest_dep + seconds2days(Saturn.T)];
[dep_SNEO,arr_SNEO,dv_SNEO,delta_t_SNEO,r1_SNEO, v1_SNEO,~,~,dv1_SNEO,dv2_SNEO] = ...
    transfer_design(dep_window_SNEO,arr_window_SNEO,Saturn.ID,NEO.ID,Sun.mu,2);
% Plots
tplot_E = linspace(dep_ES, dep_ES + seconds2days(Earth.T), 5000);
tplot_S = linspace(arr_ES, arr_ES + seconds2days(Saturn.T), 50000);
Y_E = zeros(length(tplot_E), 3); Y_S1 = zeros(length(tplot_S), 3);

for i = 1:length(tplot_E)
    [kep_E_i, ~] = uplanet(tplot_E(i), Earth.ID);
    [Y_E(i,:), ~] = kep2car(kep_E_i, Sun.mu);
end
for i = 1:length(tplot_S)
    [kep_S_i, ~] = uplanet(tplot_S(i), Saturn.ID);
    [Y_S1(i,:), ~] = kep2car(kep_S_i, Sun.mu);
end
tplot_S = linspace(dep_SNEO, dep_SNEO + seconds2days(Saturn.T), 50000);
tplot_NEO = linspace(arr_SNEO, arr_SNEO + seconds2days(NEO.T), 5000);
Y_S2 = zeros(length(tplot_S), 3); Y_NEO = zeros(length(tplot_NEO), 3);
for i = 1:length(tplot_S)
    [kep_S_i, ~] = uplanet(tplot_S(i), Saturn.ID);
    [Y_S2(i,:), ~] = kep2car(kep_S_i, Sun.mu);
end
for i = 1:length(tplot_NEO)
    [kep_NEO_i, ~, ~] = ephNEO(tplot_NEO(i), NEO.ID);
    [Y_NEO(i,:), ~] = kep2car(kep_NEO_i, Sun.mu);
end
tplot_t1 = linspace(days2seconds(dep_ES), days2seconds(arr_ES), 10000);
tplot_t2 = linspace(days2seconds(dep_SNEO), days2seconds(arr_SNEO), 10000);
y0_t1 = [r1_ES; v1_ES]; y0_t2 = [r1_SNEO; v1_SNEO];
[~, Y_t1] = ode113(@(t,y) ode_2bp(t,y,Sun.mu), tplot_t1, y0_t1, options);
[~, Y_t2] = ode113(@(t,y) ode_2bp(t,y,Sun.mu), tplot_t2, y0_t2, options);

% Earth to saturn
figure()
hold on
grid on
axis equal
plot3(Y_E(:,1), Y_E(:,2), Y_E(:,3),'Color', "#0072BD", "LineWidth", 1.5)
plot3(Y_S1(:,1), Y_S1(:,2), Y_S1(:,3),'Color', "#EDB120", "LineWidth", 1.5)
plot3(Y_t1(:,1), Y_t1(:,2), Y_t1(:,3),'Color', "#77AC30", "LineWidth", 1.5)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Earth orbit', 'Saturn orbit', 'Transfer arc','Location','best')
% Saturn to NEO
figure()
hold on
grid on
axis equal
plot3(Y_S2(:,1), Y_S2(:,2), Y_S2(:,3),'Color', "#EDB120", "LineWidth", 1.5)
plot3(Y_NEO(:,1), Y_NEO(:,2), Y_NEO(:,3),'Color', "#7E2F8E", "LineWidth", 1.5)
plot3(Y_t2(:,1), Y_t2(:,2), Y_t2(:,3),'Color', "#4DBEEE", "LineWidth", 1.5)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Saturn orbit', 'NEO orbit', 'Transfer arc','Location','best')
