%% Tracking Telemetry & Telecommand subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.25 0.25 0.3 0.3]);
addpath('..\..\commonFunctions')

%% Constants
AU = 149597870700;%m
muEarth = 3.986e14;%m3s-2
muSun = 1.327e20;%m3s-2
muJupiter = 1.267e17;%m3s-2
luminositySun = 3.846e26;%w
c = 3e8;%ms-1
magMtEarth = 7.96e15;%Tm3
magMtJupiter = 1.55e20;
Rearth = 6371;%km
Rjupiter = 69911;

%% Loading the target coordinates https://ssd.jpl.nasa.gov/horizons/app.html#/
[Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno();

%% Sources

%% Parameters
MASS = 3625;
% Define roll pitch yaw axes: X ROLL, -Y PITCH, -Z YAW
INERTIA_MOMENTS = [
    2424.8 -1.3 8.6;
    -1.3 2203.8 -0.3;
    8.6 -0.3 4521];
Imax = INERTIA_MOMENTS(3,3);
Imin = INERTIA_MOMENTS(2,2);

% TO-DO
% Find list of control modes with at least: pointing direction requirement,
% pointing accuracy or slew rate requirement and the type of control used

% Find list of Juno sensors. If possible but not strictly necessary: 
% accuracy, mass and power consumption

% Find list of Juno actuators. If possible but not strictly necessary:
% torque, mass and power consumption

%% Sizing
% Gravity gradient
Tgg.JunoToJupiter = distance.JunoToJupiter;
theta = 45; % worst case, max deviation from the local vertical
Tgg.TggEarth = 3/2*muEarth./(distance.JunoToEarth*AU).^3*(Imax - Imin)*sind(theta*2);
Tgg.TggJupiter = 3/2*muJupiter./(distance.JunoToJupiter*AU).^3*(Imax - Imin)*sind(theta*2);

% Radiation pressure
Tsrp.JunoToSun = distance.JunoToSun;
exposureArea = 1;
reflectivityIndex = 0.55;
incidence = 0; % worst case
distanceRadPressureCenters = 0.5; %Distance between the satellite's gravity and solar pressure centers
Tsrp.Tsrp = luminositySun./(4*pi*(distance.JunoToSun*AU).^2)/c*exposureArea*(1+reflectivityIndex)*cosd(incidence)*distanceRadPressureCenters;
%Albedo effect difficult to model and always smaller

% Drag
Tdrag.JunoJupiter = relVelocity.JunoJupiter;
cd = 2.2;
dragArea = 1;
distanceDragPressureCenters = 0.5;
density = 2.016e-16*exp(-(distance.JunoToEarth*AU/1000-Rearth)*10e-6/8);
Tdrag.TdragEarth = 1/2*density*cd*dragArea.*relVelocity.JunoEarth.^2;
density = 2.016e-20*exp(-(distance.JunoToJupiter*AU/1000-Rjupiter)*10e-6);
Tdrag.TdragJupiter = 1/2*density*cd*dragArea.*relVelocity.JunoJupiter.^2;

% Magnetic torque
Tmag.JunoToJupiter = distance.JunoToJupiter;
residualDipole = 10;
Tmag.TmagEarth = residualDipole*2*magMtEarth./(distance.JunoToEarth*AU).^3;
Tmag.TmagJupiter = residualDipole*2*magMtJupiter./(distance.JunoToJupiter*AU).^3;

nRCS = 4;
RCSthrust = 4;
armLength = 1.5;
maxRotation = 180;
slewRate = 0.5;
timeSlew = maxRotation/slewRate;
Fslew.x = INERTIA_MOMENTS(1,1)/nRCS/armLength*maxRotation/timeSlew^2;
Fslew.y = INERTIA_MOMENTS(2,2)/nRCS/armLength*maxRotation/timeSlew^2;
Fslew.z = INERTIA_MOMENTS(3,3)/nRCS/armLength*maxRotation/timeSlew^2

maxSlewRate = nRCS*armLength*RCSthrust/INERTIA_MOMENTS(3,3);
detumblingTimeZ = (62.5-5)*pi/180/maxSlewRate/60

mean(Tmag.TmagJupiter)*365.25*86400

%% Plots
plotOnEngine2D( ...
    time, ...
    Tgg, ...
    "Gravity gradient max contribution per body", ...
    'Distance Juno-Jupiter [AU]', ...
    'T [Nm]', ...
    'Date [year]');

plotOnEngine2D( ...
    time, ...
    Tsrp, ...
    "Solar radiation pressure", ...
    'Distance Juno-Sun [AU]', ...
    'T [Nm]', ...
    'Date [year]');

plotOnEngine2D( ...
    time, ...
    Tdrag, ...
    "Atmospheric drag", ...
    'Relative velocity Juno-Jupiter [km/s]', ...
    'T [Nm]', ...
    'Date [year]');

plotOnEngine2D( ...
    time, ...
    Tmag, ...
    "Magnetic torque", ...
    'Distance Juno-Jupiter [AU]', ...
    'T [Nm]', ...
    'Date [year]');