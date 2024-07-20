%% Electric Power subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.25 0.25 0.3 0.3]);
addpath('..\..\commonFunctions')

%% Loading the target coordinates https://ssd.jpl.nasa.gov/horizons/app.html#/
[Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno();

%% Constants
AU = 149597870700;%m
luminositySun = 3.846e26;%w

%% Juno parameters
% this page is quite good check it out https://www.lpi.usra.edu/opag/meetings/aug2015/presentations/day-2/11_beauchamp.pdf
% Juno nominal mission duration 7 years (DON'T KNOW)
duration = 10;
% EOL target power 420. Should be around 460-490 W when it arrives at Jupiter https://spaceflight101.com/juno/spacecraft-information/
powerEOL = 420;
% Working voltage 28V https://spaceflight101.com/juno/instrument-overview/
sysVoltage = 28;
% Power regulation: Fully-regulated (DON'T KNOW)
% Power control: Direct energy transfer - constant voltage (DON'T KNOW)
lineEff = 0.85; % Line efficiency for DET in daylight, from class notes
% Eclipse? Not considered. After we finish the sizing we should check the spacecraft can handle the 18 min eclipse in the fly-by though
% Class notes seem to be copied from https://ocw.mit.edu/courses/16-851-satellite-engineering-fall-2003/533672f84f6c6e3054730b672d53c87b_portfolio_nadir1.pdf

% Solar array from class notes. Find out if the solar cells are made from GaAs or Si!
arrayEffBOL = 0.316; % Efficiency BOL
arrayDeg = 0.0375; % Annual degradation of solar cell power generation
arrayflux = 100; %W/m^2 Specific power
arrayInhDeg = 0.77; % Inherent degradation of solar cell power transfer
arrayInclination = 9.7; %deg Inclination angle, assumed to be the angle between Jupiter's equatorial and Sun's ecliptic planes (IS THIS TOO CONSERVATIVE)
arraySurfaceDensity = 5.67; %kg/m^2 Surface mass density (DON'T KNOW)
%arraySurfaceDensity = 0.84
cellArea = 0.0027; %m^2 Cell area (DON'T KNOW)
cellVoltage = 2.755; %V Cell voltage (DON'T KNOW)

% Batteries
% Battery used https://www.nasa.gov/wp-content/uploads/2024/01/lithium-ion-cell-balancing-electronics-benefiting-the-satellite-industry.pdf?emrc=6645a2c301f5b
% Cobham's BEU (battery electronic unit). Each (x2) battery 3.75kg, [dm]2.92x1.33x1.33=5.17L, x8 cell, 50 Amp-hour, Lithium-ion
% The pair is Lithium-Sulphur, Li-S: https://ntrs.nasa.gov/api/citations/20210022385/downloads/NASA%20Battery_Projects_Nov2021.pdf
% Battery sizing main factor is an emergency, not the short eclipse. I've picked 3h
batteryCount = 2;
batteryTime = 3;%h
batterySpecificEnergy = 373;%W-hour/kg 746
batteryEnergyDensity = 270;%W-hour/dm^3 541
batteryCellVoltage = 3.6;%V 
DoD = 0.5; %Depth of Discharge is usually between 40 and 60%, so i picked 50%

%% Array sizing
% Solar array specific power at BOL
fluxSunlight = luminositySun./(4*pi*(max(distance.JunoToSun)*AU).^2);
fluxOutput = arrayEffBOL*fluxSunlight*arrayInhDeg;
fluxOutputBOL = fluxOutput*cosd(arrayInclination);

% Solar array specific power at EOL
lifetimeDeg = (1-arrayDeg)^duration;
fluxOutputEOL = lifetimeDeg*fluxOutputBOL;

% Power at EOL
powerRequest = powerEOL/lineEff;

% Area and mass needed
arrayArea = powerRequest/fluxOutputEOL;
arrayMass = arrayArea*arraySurfaceDensity;

% Refined sizing
Ncells = ceil(arrayArea/cellArea);
Nseries = ceil(sysVoltage/cellVoltage);
NcellsReal = ceil(Ncells/Nseries)*Nseries;

% Recalculate
arrayAreaReal = NcellsReal*cellArea;
arrayMassReal = arrayAreaReal*arraySurfaceDensity;
powerEOLreal = fluxOutputEOL*arrayAreaReal*lineEff;
fprintf("ARRAY SIZING\n")
fprintf(strcat("Array area (calculated/real): ",string(arrayAreaReal),"/","60.39"," (",string(arrayAreaReal/60.39*100-100),"%% err) m^2\n"))
fprintf(strcat("Array mass (calculated/real): ",string(arrayMassReal),"/","340"," (",string(arrayMassReal/340*100-100),"%% err) kg\n"))
fprintf(strcat("Available EOL power (calculated/real): ",string(powerEOLreal),"/","420"," (",string(powerEOLreal/420*100-100),"%% err) W\n\n"))

%% Power demanded
ellapsed.year = time.year-time.year(1);
availablePower.SunDistance = distance.JunoToSun;

fluxSunlight = luminositySun./(4*pi*(distance.JunoToSun*AU).^2);
fluxOutput = arrayEffBOL*fluxSunlight*arrayInhDeg;
fluxOutputBOL = fluxOutput*cosd(arrayInclination);
lifetimeDeg = (1-arrayDeg).^ellapsed.year;
fluxOutputReal = lifetimeDeg.*fluxOutputBOL;
availablePower.Available = fluxOutputReal*arrayAreaReal*lineEff;
availablePower.Available(availablePower.Available > 800) = 800;

Q_int = ones(size(distance.JunoToSun));
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


% TCS params
sigma=5.67e-8; %Stefan-Boltzman constant [W/m^2/K^4]
T_ds=3; %deep space temperature

eqSphereR = sqrt(7.6/pi);
AcrossSection = pi*eqSphereR^2;
Atotal = 4*pi*(eqSphereR)^2;
alpha_sc=0.05; %absorptivity of thermal coating
eps_sc=0.05; %emissivity of thermal coating
Tmin = 273 - 5;
%
Q_sun = AcrossSection*fluxSunlight*alpha_sc;
Q_emitted = Q_int+Q_sun;
Q_heaters = sigma*eps_sc*Atotal*(Tmin^4-T_ds^4) - Q_emitted;
Q_heaters(Q_heaters < 0) = 0;


availablePower.Demanded = Q_int+Q_heaters;

%% Battery sizing Sizing is ok but the numbers are wildly off here
batteryCapacity = batteryTime*max(availablePower.Demanded)/DoD/batteryCount/lineEff;
batteryMass = batteryCapacity/batterySpecificEnergy;
batteryVolume = batteryCapacity/batteryEnergyDensity;

% Refined sizing
NcellsBattery = ceil(sysVoltage/batteryCellVoltage);
voltageReal = NcellsBattery*batteryCellVoltage;

batteryPackageEff = 0.8;
batteryCellCapacity = 10;
stringCapacity = batteryPackageEff*batteryCellCapacity*voltageReal;

NparallelStrings = ceil(batteryCapacity/stringCapacity);
batteryCapacityReal = NparallelStrings*stringCapacity;
batteryCapacityReal = batteryCapacityReal/voltageReal; %amp hour

fprintf("BATTERY SIZING (there are 2 batteries)\n")
fprintf(strcat("Battery cells (calculated/real): ",string(NcellsBattery),"/","8"," (",string(NcellsBattery/8*100-100),"%% err)\n"))
fprintf(strcat("Battery voltage (calculated/real): ",string(voltageReal),"/","28"," (",string(voltageReal/28*100-100),"%% err) V\n"))
fprintf(strcat("Battery capacity (calculated/real): ",string(batteryCapacityReal),"/","50"," (",string(batteryCapacityReal/50*100-100),"%% err) Amp-hour\n"))
fprintf(strcat("Battery mass (calculated/real): ",string(batteryMass),"/","3.75"," (",string(batteryMass/3.75*100-100),"%% err) kg\n"))
fprintf(strcat("Battery volume (calculated/real): ",string(batteryVolume),"/","5.17"," (",string(batteryVolume/5.17*100-100),"%% err) L\n"))

%% Plots
plotOnEngine2D( ...
    ellapsed, ...
    availablePower, ...
    'Available vs. demanded power', ...
    'Distance Juno-Sun [AU]', ...
    'Internal power [W]', ...
    'Ellapsed time since BOL [years]');