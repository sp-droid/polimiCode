%% Propulsion subsystem sizing
clc

%% Key
% X === Accurate value found
% R === Accurate range found, but handpicked a value *using other missions*
% G === Complete guess, *using other missions*
% Nothing === Pending

%% Parameters
%MAIN_ENGINE_DELTA_V = 1646.9;           %X [m/s] Delta V budget of the biprop maneuvers IMPORTANT
%RCS_DELTA_V = 385.5;                    % ídem for rcs  IMPORTANT
%X Delta-V budget breakdown by maneuver and engine
deltaVbreakdown = [4.4, 754.9, 31.1, 892, 350];
deltaVengineType = [ 1,     0,    1,   0,   1];
MAIN_ENGINE_ISP = 317;                  %X [s] Isp of the main engine
RCS_ISP = 229;                          %R ídem for rcs
M_DRY = 1593;                           %X [kg] Dry mass at launch

OF_RATIO = 0.85;                        %R [-] OxFu ratio selected
RHO_FU = 1008;                          %X [kg/m3] Fuel density
RHO_OX = 1447;                          %X [kg/m3] Oxidizer density

P_INLET = 21e5;                         %R [Pa] Combustion chamber inlet pressure
P_INJ_DROP = 0.3;                       % [-] Injection pressure drop as fraction of P_CC
P_FEED_DROP = 50000;                    % [Pa] Feeding lines pressure drop

GAMMA_PRESS = 1.67;                     %X [-] Pressurant gamma value
R_GAS_PRESS = 2077.3;                   %X [J/kgK] Pressurant specific gas constant
LIQUID_DENSITY_PRESS = 125;             %X [kg/m3] Liquid helium density
P_PRESS_F = 1;                          %G [-] ? pressure as fraction of pTank
P_PRESS_I = 10;                         %G [-] Pressurant tank pressure as fraction of pTank
T_TANKS = 293;                          %G [K] Temperature of stored propellant
%T_PRESS = 293;                         % [K] Temperature of stored pressurant Not needed anymore since we assume liquid helium

N_TANKS_FU = 4;                         %X [-] Number of fuel tanks
N_TANKS_OX = 2;                         %X [-] Number of oxidizer tanks
N_TANKS_PRESS = 2;                      %X [-] Number of pressurant tanks
% Shape assumed to be spherical
TANK_PROP_RHO = 4430;                   %G [kg/m3] Density of material used for the tank wall (Ti4Al6Vi)
TANK_PROP_SIGMA = 880e6;                %G [Pa] Uniaxial yield stress of (idem)
TANK_PRESS_RHO = 2700;                  %G [kg/m3] Density of material used for the pressurant tank wall 
TANK_PRESS_SIGMA = 285e6;               %G [Pa] Uniaxial yield stress of (idem)
TANK_PRESS_LR_RATIO = 4;                %G [-] Length/Radius ratio for the pressurant cylinder tanks

MASS_MAIN_ENGINE = 4.5;                 %X [kg] Mass of main engine
MASS_RCS_THRUSTER = 0.33;               %X [kg] Mass of thrusters
N_RCS_THRUSTER = 12;                    %X [-] Number of thrusters

POWER_MAIN_VALVES = 0;                  % For power budget calculations
POWER_PER_RCS_VALVE = 13.64;            %X [W] For (max) power budget calculations
POWER_FEEDLINES = 0;                    % For power budget calculations

%% Real targets
REAL_PROPELLANT_MASS = 2032;            %X [kg]
REAL_PS_MASS = 0;                       % [kg]

%% Constants
G0 = 9.80665;

%% Sizing
mProp = 0; mPropRCS = 0;
for i=length(deltaVbreakdown):-1:1
    if deltaVengineType(i)==0
        mProp = mProp + ((1.20)*M_DRY + mProp + mPropRCS) * (exp(1.00*deltaVbreakdown(i) / MAIN_ENGINE_ISP / G0) - 1);
    else
        mPropRCS = mPropRCS + ((1.20)*M_DRY + mProp + mPropRCS) * (exp(1.00*deltaVbreakdown(i) / RCS_ISP / G0) - 1);
    end
end

mProp = mProp*(1.02+0.03+0.005); mPropRCS = mPropRCS*(1.02+0.03+0.005);
% Propellant mass, deltaV = Isp*g0*ln(m0/mFinal)
%mProp = (1.02+0.03+0.005)*(1.2)*M_DRY * (exp((1.1)*MAIN_ENGINE_DELTA_V / MAIN_ENGINE_ISP / G0) - 1);
%mPropRCS = (1.02+0.03+0.005)*(1.2)*M_DRY * (exp((1.1)*RCS_DELTA_V / RCS_ISP / G0) - 1);
fprintf(strcat("Propellant for main engine: ",string(mProp)," kg\n"));
fprintf(strcat("Extra fuel for RCS: ",string(mPropRCS)," kg\n"));
fprintf(strcat("Total: ",string(mProp+mPropRCS)," kg;  Real: ",string(REAL_PROPELLANT_MASS)," kg;  Difference: ",string((mProp+mPropRCS)/REAL_PROPELLANT_MASS*100-100),"%%\n"));

% OxFu masses (assuming equal tank volumes, OF = mOx / mFu = rhoOx / rhoFu
mFu = mProp / (1 + OF_RATIO);
mOx = mProp - mFu;
mFu = mFu + mPropRCS;
volFu = (1.1)*mFu / RHO_FU;
volOx = (1.1)*mOx / RHO_OX;
volProp = volFu + volOx;
fprintf(strcat("Fuel mass: ",string(mFu)," kg;\t Oxidizer mass: ",string(mOx)," kg\n"));
fprintf(strcat("Fuel volume: ",string(volFu*1000)," L;\t Oxidizer volume: ",string(volOx*1000)," L\n"));

% Fuel/Oxidizer tank pressure
pTank = P_INLET + P_FEED_DROP;
pCC = P_INLET / (1+P_INJ_DROP);
fprintf(strcat("Combustion chamber pressure: ",string(pCC/100000)," bar\n"));
fprintf(strcat("Combustion chamber inlet pressure: ",string(P_INLET/100000)," bar\n"));
fprintf(strcat("Propellant tank pressure: ",string(pTank/100000)," bar\n"));
% Pressurant tank pressure, mass and volume
pPressF = P_PRESS_F * pTank;
pPressI = P_PRESS_I * pTank;
fprintf(strcat("Pressurant tank pressure: ",string(pPressI/100000)," bar\n"));
mPress = pTank * volProp / R_GAS_PRESS / T_TANKS * GAMMA_PRESS / (1 - pPressF/pPressI);
fprintf(strcat("Pressurant mass: ",string(mPress)," kg\n"));
%volPress = mPress * R_GAS_PRESS * T_PRESS / pPressI
volPress = mPress / LIQUID_DENSITY_PRESS;
fprintf(strcat("Pressurant volume: ",string(volPress*1000)," L\n"));

% Tanks radius, thickness and empty mass
rFu = (3/4/pi*volFu/N_TANKS_FU)^(1/3);
rOx = (3/4/pi*volOx/N_TANKS_OX)^(1/3);
rPress = (1/TANK_PRESS_LR_RATIO/pi*volPress/N_TANKS_PRESS)^(1/3);
lPress = TANK_PRESS_LR_RATIO*rPress;
tFu = pTank * rFu / 2 / TANK_PROP_SIGMA;
tOx = pTank * rOx / 2 / TANK_PROP_SIGMA;
tPress = pPressI * rPress / TANK_PRESS_SIGMA;
mTankFu = TANK_PROP_RHO * 4/3*pi*((rFu+tFu)^3 - rFu^3);
mTankOx = TANK_PROP_RHO * 4/3*pi*((rOx+tOx)^3 - rOx^3);
mTankPress = TANK_PRESS_RHO * lPress*pi*((rPress+tPress)^2 - rPress^2);
fprintf(strcat("Fuel       tank; Radius: ",string(rFu)," m;  Thickness: ",string(tFu*1000)," mm;  Mass: ",string(mTankFu)," kg; Quantity: ",string(N_TANKS_FU),"\n"));
fprintf(strcat("Oxidizer   tank; Radius: ",string(rOx)," m;  Thickness: ",string(tOx*1000)," mm;  Mass: ",string(mTankOx)," kg; Quantity: ",string(N_TANKS_OX),"\n"));
fprintf(strcat("Pressurant tank; Radius: ",string(rPress)," m; Length: ",string(lPress)," m;  Thickness: ",string(tPress*1000)," mm;  Mass: ",string(mTankPress)," kg; Quantity: ",string(N_TANKS_PRESS),"\n"));

massPS = (1.1) * (mTankFu*N_TANKS_FU + mTankOx*N_TANKS_OX + mTankPress*N_TANKS_PRESS + mPress + MASS_RCS_THRUSTER*N_RCS_THRUSTER + MASS_MAIN_ENGINE);
fprintf(strcat("Propulsion system mass budget: ",string(massPS)," kg\n"));

% Power budget
powerPS = POWER_MAIN_VALVES + POWER_PER_RCS_VALVE*N_RCS_THRUSTER + POWER_FEEDLINES;
fprintf(strcat("Propulsion system power budget: ",string(powerPS)," W\n"));