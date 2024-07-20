%% Thermal Control subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.25 0.25 0.3 0.3]);
addpath('..\..\commonFunctions')

%% Constants
AU = 149597870700;%m

%% Loading the target coordinates https://ssd.jpl.nasa.gov/horizons/app.html#/
[Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno();

%% Ray intersection (eclipse conditions for Sun-Earth-Juno and Sun-Jupiter-Juno systems)
eclipse.JunoEarth = eclipseCondition(Earth, Juno, 6371*1000/AU);
eclipse.JunoJupiter = eclipseCondition(Jupiter, Juno, 69911*1000/AU);
% figure;
% scatter(time.year, eclipse.JunoEarth)
% figure;
% scatter(time.year, eclipse.JunoJupiter)

%% Constants
luminositySun = 3.846e26;%w
sigma=5.670367e-8; % Stefan-Boltzman constant [W/m^2/K^4]
T_ds=3; % deep space temperature

% SHAPE OF JUNO
% Juno is an hexagonal prism with the HGA on top and 3 big solar panels
%Atop = 7.6;

eqSphereR = sqrt(7.6/pi);
fprintf(strcat("Equivalent sphere radius: ",string(eqSphereR)," m\n"))
AcrossSection = pi*eqSphereR^2;
Atotal = 4*pi*(eqSphereR)^2;
fprintf(strcat("Sphere cross section/total area: ",string(AcrossSection)," / ",string(Atotal)," m^2\n"))
AsolarPanels = 59;

Q_int_hot = 171.78;
Q_int_cold = 259.97; % higher, in the cold case (at Jupiter) more instruments working at the same time

a_jupiter = 0.41;
T_Jupiter = 128;
eps_Jupiter = 0.27; % https://www.sciencedirect.com/science/article/pii/0019103565900400#:~:text=A%20simple%20interpretation%20of%20the,K%20and%20emissivity%20of%200.27.
R_Jupiter = 69911;

% Parameters
alpha_sc=0.04; % absorptivity of thermal coating
eps_sc=0.04; % emissivity of thermal coating
alpha_panel = 0.7;
eps_panel = 0.7;

eps_radiator = 0.74;


%% SIZING - WORST HOT CASE: Deep space 0.8 AU

fprintf(strcat("--------HOT CASE-------- Q_int: ",string(Q_int_hot)," W\n"))
q_sun = luminositySun./(4*pi*(min(distance.JunoToSun)*AU).^2);
% main body
Q_sun = AcrossSection*q_sun*alpha_sc;
fprintf(strcat("qSun / Qsun: ",string(q_sun)," / ",string(Q_sun)," W\n"))
Q_emitted = Q_int_hot+Q_sun;
T_sc = (Q_emitted/(sigma*eps_sc*Atotal)+T_ds^4)^(1/4);
fprintf(strcat("Spacecraft temperature hot case: ",string(T_sc-273)," ºC\n"))

Tmax = 273 + 30;
Q_radiator = Q_emitted - sigma*eps_sc*Atotal*(Tmax^4-T_ds^4);
fprintf(strcat("Spacecraft radiator power for hot case and Tmax=",string(Tmax-273)," ºC: ",string(Q_radiator)," W \n"))
Aradiator = Q_radiator/(sigma*eps_radiator)/(Tmax^4-T_ds^4);
fprintf(strcat("Spacecraft radiator area for hot case and Tmax=",string(Tmax-273)," ºC: ",string(Aradiator)," m^2 \n"))

% panels
Q_sun = AsolarPanels*q_sun*alpha_panel;
fprintf(strcat("qSun / Qsun: ",string(q_sun)," / ",string(Q_sun)," W\n"))
Q_emitted = Q_sun-Q_int_hot;
T_panels = (Q_emitted/(sigma*eps_panel*AsolarPanels*2)+T_ds^4)^(1/4);
fprintf(strcat("Panel temperature hot case: ",string(T_panels-273)," ºC\n"))

%% SIZING - COLD CASE: Eclipse behind Jupiter
fprintf(strcat("--------COLD CASE-------- Q_int: ",string(Q_int_cold)," W\n"))
q_sun = luminositySun./(4*pi*(max(distance.JunoToSun)*AU).^2);
q_sun_sc = 1367.5*(1/max(distance.JunoToSun))^2;
q_albedo = q_sun*a_jupiter*(R_Jupiter/(74111))^2;
T_Jupiter = ((1-a_jupiter)*q_sun/(sigma*eps_Jupiter*4))^(1/4);
q_IR = eps_Jupiter*sigma*T_Jupiter^4*(R_Jupiter/74111)^2;

%main body
Q_sun = AcrossSection*q_sun*alpha_sc;
fprintf(strcat("qSun / Qsun: ",string(q_sun)," / ",string(Q_sun)," W\n"))
Q_albedo = AcrossSection*q_albedo*alpha_sc;
Q_IR = AcrossSection*q_IR*alpha_sc;
Q_emitted = Q_int_cold+Q_sun+Q_albedo+Q_IR;
T_sc = (Q_emitted/(sigma*eps_sc*Atotal)+T_ds^4)^(1/4);
fprintf(strcat("Spacecraft temperature cold case: ",string(T_sc-273)," ºC\n"))

Tmin = 273 + 10;
Q_heaters = sigma*eps_sc*Atotal*(Tmin^4-T_ds^4) - Q_emitted;
fprintf(strcat("Spacecraft heater power for cold case and Tmin=",string(Tmin-273)," ºC: ",string(Q_heaters)," W \n"))

% panels - Thermally decoupled from the rest of the spacecraft?
Q_sun = AsolarPanels*q_sun*alpha_panel;
Q_albedo = AsolarPanels*q_albedo*alpha_panel;
Q_IR = AsolarPanels*q_IR*alpha_panel;

fprintf(strcat("qSun / Qsun: ",string(q_sun)," / ",string(Q_sun)," W\n"))
Q_emitted = Q_sun+Q_IR+Q_albedo-Q_int_cold;
T_panels = (Q_emitted/(sigma*eps_panel*AsolarPanels*2)+T_ds^4)^(1/4);
fprintf(strcat("Panel temperature cold case: ",string(T_panels-273)," ºC\n"))
Q_emitted = Q_emitted - Q_heaters;
T_panels = (Q_emitted/(sigma*eps_panel*AsolarPanels*2)+T_ds^4)^(1/4);
fprintf(strcat("Panel temperature cold case + accounting for heaters: ",string(T_panels-273)," ºC\n"))

fprintf(strcat("Max power allowed: 435 W. Internal + heaters: ",string(Q_heaters+Q_int_cold)," W\n"))

%% GRAPH
q_sun = luminositySun./(4*pi*(distance.JunoToSun*AU).^2);
Q_sun = AcrossSection*q_sun*alpha_sc;
Q_int = Q_int_hot*ones(size(distance.JunoToSun));
for i=1:length(Q_int)
    if (time.year(i) < 2013.332)
        Q_int(i) = 186.56;
    elseif (time.year(i) < 2013.76)
        Q_int(i) = 171.78;
    elseif (time.year(i) < 2013.77)
        Q_int(i) = 260.96;
    elseif (time.year(i) < 2016.33)
        Q_int(i) = 171.78;
    elseif (time.year(i) < 2016.80)
        Q_int(i) = 259.97;
    elseif (time.year(i) < 2016.95)
        Q_int(i) = 304.76;
    elseif (time.year(i) < 2017.09)
        Q_int(i) = 259.97;
    elseif (time.year(i) < 2017.67)
        Q_int(i) = 304.76;
    else
        Q_int(i) = 259.97;
    end
end
Q_emitted = Q_int+Q_sun;
T_sc = (Q_emitted./(sigma*eps_sc*Atotal)+T_ds^4).^(1/4);
tcs.Temperature = T_sc;
tcs.Temperature(tcs.Temperature > Tmax) = Tmax;
tcs.Temperature(tcs.Temperature < Tmin) = Tmin;
tcs.Temperature = tcs.Temperature-273;

Q_radiator = Q_emitted - sigma*eps_sc*Atotal*(Tmax^4-T_ds^4);
Q_radiator(Q_radiator < 0) = 0;
Aradiator = Q_radiator/(sigma*eps_radiator)/(Tmax^4-T_ds^4);
Aradiator = Aradiator/(0.16*3)*100;
tcs.Radiators = Aradiator;

Q_heaters = sigma*eps_sc*Atotal*(Tmin^4-T_ds^4) - Q_emitted;
Q_heaters(Q_heaters < 0) = 0;
Q_heaters = Q_heaters./(435-Q_int)*100;
tcs.Heaters = Q_heaters;

plotOnEngine2D( ...
    time, ...
    tcs, ...
    "TCS percentage usage (100 % = maximum allowed)", ...
    'Spacecraft main body temperature [C]', ...
    'Usage [%]', ...
    'Date [year]');

qint.JunoToJupiter = distance.JunoToJupiter;
qint.Qint = Q_int;
plotOnEngine2D( ...
    time, ...
    qint, ...
    "Internal power", ...
    'Distance Juno-Jupiter [AU]', ...
    'Power [W]', ...
    'Date [year]');
function condition = eclipseCondition(planetPos, satellitePos, planetRadius)

R2 = planetRadius^2;
r = vecnorm(planetPos(:,1:3), 2, 2); % Distance Sun-Planet
s = vecnorm(satellitePos(:,1:3), 2, 2); % Distance Sun-Satellite
tc = dot(planetPos(:,1:3),satellitePos(:,1:3), 2)./s; % Projection of the planet position vector into the satellite's
d2 = r.^2 - tc.^2;

condition = (d2 < R2) & (s > r) & (tc > 0);
end