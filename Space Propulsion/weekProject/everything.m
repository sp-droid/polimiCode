clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Constants
g0 = 9.80665; %[m/s^2]
R = 8.314; % [J/mol/K]
T0 = 298; %[K]

%% Given parameters
deltaV = 2500;
massFinal = 250;

%% Case
propellants = 'green';
scaling = 'doubled';
out = highleveldecisions(propellants, scaling, R, g0);

%% General calculations
% Total propellant mass [kg]
massProp = massFinal*exp(deltaV/out.Isp/g0)-massFinal
% Propellants mass flow rates [kg/s]
mFU = out.mP/(out.OF+1)
mOX = out.mP-mFU
% Total burn time [min]
ToF = massProp/out.mP/60

%% Geometry
% Main diameters
D = sqrt(4*out.A/pi);

% CC volume
Vc = out.Lstar*out.A(3);
% Completing the cylinder
Lcc = 1/out.contraction*(Vc/out.A(3)-1/3*sqrt(out.A(3)/pi)/tand(out.convergentAngle)*(out.contraction^(1/3)-1));

% Conical nozzle
Lcon = 1/tand(out.convergentAngle)*(D(1)-D(3))/2;
Ldiv = 1/tand(out.divergentAngle)*(D(4)-D(3))/2;
Ltotal = Lcc+Lcon+Ldiv;

% PLOT geometry
n = 1000;
xi = linspace(0,Lcc+Lcon+Ldiv,n)';
deltax = xi(2)-xi(1);
Di = zeros(n,1); Ai = zeros(n,1);
for i=1:n
    Di(i) = diamEngineConical(xi(i), Lcc, Lcon, D(1), D(3), out.convergentAngle, out.divergentAngle);
    Ai(i) = pi*Di(i)^2/4;
end

figure();
plot(xi*100,Di/2*100,'LineWidth',2,'Color','black')
hold on
plot(xi*100,-Di/2*100,'LineWidth',2,'Color','black')
hold on
plot([0,0],[-D(1),D(1)]/2*100,'LineWidth',2,'Color','black','LineStyle','-')
hold on
plot([Ltotal,Ltotal]*100,[-D(4),D(4)]/2*100,'LineWidth',2,'Color','black','LineStyle','--')
hold on
plot([0,Ltotal]*100,[0,0],'LineWidth',2,'Color','black','LineStyle','-.')
xlabel('\boldmath$L\hspace{0.5em}[cm]$','Interpreter','latex'); ylabel('\boldmath$Y\hspace{0.5em}[cm]$','Interpreter','latex');
title('\boldmath$Cross-section$','Interpreter','latex')
set(gca,'fontsize', 16) 
axis equal;
grid on;
hold off

graphy.A = Ai/out.A(3);
plotOnEngine2D(graphy,'\boldmath$Area$','\boldmath$A/A_{t}\hspace{0.5em}[-]$',xi,Di,D,Ltotal); clearvars graphy

%% Solve for each element
% AreaPrandtl, Reynolds
Mi = zeros(n,1); Pi = zeros(n,1); Ti = zeros(n,1); rhoi = zeros(n,1); ai = zeros(n,1); vi = zeros(n,1);
mui = zeros(n,1); Rei = zeros(n,1); Pri = zeros(n,1);
for i=1:n
    res = solveEngine(Ai(i), xi(i), Lcc+Lcon, out); %[M,P,T,rho,a,v]
    Mi(i) = res(1); Pi(i) = res(2); Ti(i) = res(3); rhoi(i) = res(4); ai(i) = res(5); vi(i) = res(6);
    % Sutherland's law https://es.wikipedia.org/wiki/Ley_de_Sutherland (replace wiki reference!)
    C = 120;
    mui(i) = out.mu*(Ti(i)/out.T(1))^1.5*(out.T(1)+C)/(Ti(i)+C);
    Rei(i) = rhoi(i)*vi(i)*Di(i)/mui(i);
    Pri(i) = out.cp*mui(i)/out.K;
end
% PLOT properties
graphy.M = Mi;
plotOnEngine2D(graphy,'\boldmath$Mach\hspace{0.5em}number$','\boldmath$M\hspace{0.5em}[-]$',xi,Di,D,Ltotal); clearvars graphy

%% Heat transfer
Tawi = zeros(n,1); Nui = zeros(n,1); qGi = zeros(n,1);
for i=1:n
    Ttotal = Ti(i)*(1+0.5*(out.gamma-1)*Mi(i)^2);
    RecoveryFactor = (1+Pri(i)^(1/3)*0.5*(out.gamma-1)*Mi(i)^2)/(1+0.5*(out.gamma-1)*Mi(i)^2);
    Tawi(i) = RecoveryFactor*Ttotal;

    % Dittus Boelter. Twall is constrained by the maximum temperature of the material
    %Nui(i) = 0.0265*Rei(i)^(4/5)*Pri(i)^(0.3);
    %hG = Nui(i)*out.K/Di(i);

    % Bartz model https://ntrs.nasa.gov/api/citations/19990025912/downloads/19990025912.pdf
    sigma = 1/(0.5*out.material.T/Tawi(i)*(1+(out.gamma-1)/2*Mi(i)^2)+0.5)^0.68/(1+(out.gamma-1)/2*Mi(i)^2)^0.12;
    Rcurv = 0.012;
    hG = sigma*(out.A(3)/Ai(i))^0.9*(0.026/D(3)^0.2*mui(i)^0.2*out.cp/Pri(i)^0.6*(out.P(1)/out.cstar)^0.8*(D(3)/Rcurv)^0.1);

    hG = hG/(1+out.sootEffect.*hG);
    qGi(i) = hG*(Tawi(i)-out.material.T);
end
graphy.q = qGi;
plotOnEngine2D(graphy,'\boldmath$Heat\hspace{0.5em}flux$','\boldmath$\dot{q}\hspace{0.5em}[\frac{W}{m^{2}s}]$',xi,Di,D,Ltotal); clearvars graphy

%% CC heat transfer
Qcooling = pi*Di(1)*deltax*qGi(1); qRingi = zeros(n,1);
for i=2:n
    Acooling = pi*(Di(i)+Di(i-1))/2*deltax;
    Qcooling = Qcooling + Acooling*(qGi(i)+qGi(i-1))/2;
    % Instead of the heat flux in each point, this is on each section
    qRingi(i) = pi*Di(i)*qGi(i);
end
graphy.qRing = qRingi;
plotOnEngine2D(graphy,['\boldmath$Section\hspace{0.3em}heat\hspace{0.3em}flux\hspace{0.3em}(total=',sprintf('%.2f',Qcooling/1000),'kW)$'],'\boldmath$\dot{q}\hspace{0.5em}[\frac{W}{ms}]$',xi,Di,D,Ltotal); clearvars graphy
Qcooling

%% Cooling
Pregen = 11e5;
% Possible regenerative cooling using fuel
Tliquid = out.lOXboilT(Pregen)
deltaTox = Qcooling/mOX/out.lOXcp(Tliquid)
ToxBudget = Tliquid-deltaTox-T0
QoxMax = (Tliquid-T0)*mOX*out.lOXcp(Tliquid)
% Possible regenerative cooling using oxidizer
Tliquid = out.lFUboilT(Pregen)
deltaTfu = Qcooling/mFU/out.lFUcp(Tliquid)
TfuBudget = Tliquid-deltaTfu-T0
QfuMax = (Tliquid-T0)*mFU*out.lFUcp(Tliquid)

Qregen = max(QoxMax,QfuMax);

mfilm = 0; mfilmi = zeros(n,1);
roughnessHeight = 25e-6;
filmEfficiency = 0.7;
stopRegen = false;
options = optimset('Display','off');
for i=1:n
    Acooling = pi*Di(i)*deltax;
    Qregen = Qregen - Acooling*qGi(i);
    if Qregen>0
        continue
    end
    if stopRegen==false
        disp('Regen vs film point')
        disp(xi(i))
        disp(Lcc)
        disp(Lcc+Lcon)
        disp(Ltotal)
    end
    stopRegen = true;

    % Colebrook-White equation to get Darcy's coefficient
    func = @(x) cWDarcyCoeff(x, Di(i), Rei(i), roughnessHeight);
    f = fsolve(func, 0.1, options);
    
    a = 2*1/f;
    cpgas = out.gFUcp(out.material.T);
    cpliq = out.lFUcp((out.lFUboilT(Pi(i))+T0)/2);
    H = cpgas*(Tawi(i)-out.material.T)/(cpliq*(out.material.T-T0)+out.FUvapEnthalpy);
    mfilmi(i) = (out.mP+mfilm)*H/filmEfficiency/a/(1+0);
    mfilm = mfilm + Acooling*mfilmi(i);

%     if abs(xi(i)-Ltotal+0.0001)<=0.0005
%         disp(Tawi(i))
%         disp(out.material.T)
%         disp(f)
%         disp(cpgas)
%         disp(cpliq)
%         disp(mfilmi(i))
%     end
end
graphy.mfilm = mfilmi*1000;
plotOnEngine2D(graphy,['\boldmath$Film\hspace{0.3em}cooling\hspace{0.3em}pointwise\hspace{0.3em}mass\hspace{0.3em}flow-rate\hspace{0.3em}(total=',sprintf('%.2f',mfilm*1000),'g/s)$'],'\boldmath$\dot{m}_{film}\hspace{0.5em}[\frac{g}{m^{2}s}$]',xi,Di,D,Ltotal); clearvars graphy
mfilm

%% Stresses
tWi = zeros(n,1); tWmaxregenOXi = zeros(n,1); tWmaxregenFUi = zeros(n,1);
for i=1:n
    % Temporary measure of stress in cylinder
    tWi(i) = Pi(i)*Di(i)/2/out.material.yieldAtT;
    tWmaxregenOXi(i) = out.material.k*(out.material.T-out.lOXboilT(Pregen))/qGi(i);
    tWmaxregenFUi(i) = out.material.k*(out.material.T-out.lFUboilT(Pregen))/qGi(i);
end
%graphy.maxThicknessRegenOX = tWmaxregenOXi*1000;
%graphy.maxThicknessRegenFU = tWmaxregenFUi*1000;
%plotOnEngine2D(graphy,'Maximum thickness','tw [mm]',xi,Di,D,Ltotal); clearvars graphy
graphy.PressureStress = tWi*1000;
plotOnEngine2D(graphy,'\boldmath$Minimum\hspace{0.5em}wall\hspace{0.5em}thickness$','\boldmath$t_{w}\hspace{0.5em}[mm]$',xi,Di,D,Ltotal); clearvars graphy

%% Injectors
AinjFU = pi*out.injFUd^2/4;
AinjOX = pi*out.injOXd^2/4;
ratioAreaInj = (AinjFU*out.injFUn+AinjOX*out.injOXn)/(pi*D(1)^2/4)
ratioPinjFU = (mFU/out.injFUn/out.injCd/AinjFU)^2/2/out.lFUrho/out.P(1)
ratioPinjOX = (mOX/out.injOXn/out.injCd/AinjOX)^2/2/out.lOXrho/out.P(1)
VinjFU = out.injCd*sqrt(2*ratioPinjFU*out.P(1)/out.lFUrho)
VinjOX = out.injCd*sqrt(2*ratioPinjOX*out.P(1)/out.lOXrho)
mFU
mOX

if strcmp(propellants,'toxic')
    % Like-on-like impinging injectors, doublets
    % The angle can be the same for both jets in each doublet, 45º?
elseif strcmp(propellants,'green')
    % Unlike impinging injectors. Pentads
    % Angle needs adjustment based on jet momentum difference
end

%% Functions
function out = highleveldecisions(propellants, scaling, R, g0)
% eps, Ae/At [-]

% contr, contraction ratio Ac/At [-]
% Lstar, characteristic chamber length [m]

% OF, ox/fuel ratio [-]
% P, pressure [Pa]
% T, temperature [K]
% rho, density [kg/m^3]
% Mm, molecular mass [kg/mol]
% gamma, heat capacity ratio [-]
% cstar, characteristic velocity [m/s]
% Isp, specific impulse [m/s]
% cp, specific heat [J/kg/K], frozen
% mu, viscosity [Pa*s]
% k, conductivity [W/m/K]

% F, thrust [newtons]

% Lstar reference (swap for scientific one): https://space.stackexchange.com/questions/31009/how-are-the-combustion-chamber-length-and-diameter-decided

% In vectors, first is CC, second is throat, third is nozzle exit
out.eps = 200; % DECISION
out.convergentAngle = 45; % DECISION [deg]
out.divergentAngle = 15; % DECISION [deg]
% P = 10; % DECISION, included in the CEA results
if strcmp(propellants,'toxic')
    out.contraction = 5.90; % DECISION
    out.Lstar = 0.8255; % DECISION
    out.sootEffect = 0;

    out.material.name = 'Ti-6Al-4V'; % DECISION
    out.material.yieldAtT = 580*10^6;
    out.material.T = 700;
    out.material.k = 6.7;
    
    %https://www.chemeurope.com/en/encyclopedia/
    %https://webbook.nist.gov/cgi/
    %https://web.stanford.edu/group/haiwanglab/HyChem/approach/Report_RP2_Fuel_Thermochemical_Properties_v2.pdf
    %http://imartinez.etsiae.upm.es/~isidoro/dat1/eLIQ.pdf
    out.lFUmm = 46.07/1000;
    out.lFUcp = @(T) cpMMH(T);
    out.lFUboilT = @(P) TboilMMH(P);
    out.lFUrho = 874;
    out.lFUmu = 0.97/1000;
    out.lOXmm = 92.011/1000;
    out.lOXcp = @(T) cpN2O4(T);
    out.lOXboilT = @(P) TboilN2O4(P);
    out.lOXrho = 1431;
    out.lOXmu = 0.47/1000;

    out.FUvapEnthalpy = 8.7621e+05;
    out.gFUcp = @(T) cpGasMMH(T);

    out.P = [1000000;0;0;0];
    out.OF = 1.69;
    out.T = [3079.67;0;0;0];
    out.rho = [0.80499;0;0;0];
    out.Mm = 20.613/1000;
    out.M = [0.000;0.101;1.000;0];
    
    out.cp = 2.1472*1000;
    out.mu = 0.96015/10000;
    out.K = 3.2823/10;

    out.injCd = 0.7;
    out.injFUd = 0.001;
    out.injOXd = 0.001;
    if strcmp(scaling,'halved')
        out.F = 480/2; % DECISION
        out.injFUn = 3; % DECISION
        out.injOXn = 4; % DECISION
    elseif strcmp(scaling,'nominal')
        out.F = 480; % DECISION
        out.injFUn = 6; % DECISION
        out.injOXn = 8; % DECISION
    elseif strcmp(scaling,'doubled')
        out.F = 480*2; % DECISION
        out.injFUn = 12; % DECISION
        out.injOXn = 16; % DECISION
    end
elseif strcmp(propellants,'green')
    out.contraction = 5.93; % DECISION
    out.Lstar = 1.651; % DECISION
    out.sootEffect = 0;

    out.material.name = 'INCONEL X-750'; % DECISION
    out.material.yieldAtT = 700*10^6;
    out.material.T = 922;
    out.material.k = 20.6;
    
    out.lFUmm = 167.9/1000;
    out.lFUcp = @(T) cpRP1(T);
    out.lFUboilT = @(P) TboilRP1(P);
    out.lFUrho = 820;
    out.lFUmu = 2.4/1000;
    out.lOXmm = 34.014/1000;
    out.lOXcp = @(T) cpH2O2(T);
    out.lOXboilT = @(P) TboilH2O2(P);
    out.lOXrho = 1450;
    out.lOXmu = 1.2/1000;
    
    out.FUvapEnthalpy = 2.2635e+05;
    out.gFUcp = @(T) cpGasRP1(T);

    out.P = [1000000;0;0;0];
    out.OF = 5.88;
    out.T = [2810.88;0;0;0];
    out.rho = [0.90809;0;0;0];
    out.Mm = 21.223/1000;
    out.M = [0.000;0.101;1.000;0];

    out.cp = 2.5031*1000; %frozen
    out.mu = 0.98040/10000;
    out.K = 3.3851/10;

    out.injCd = 0.7;
    out.injFUd = 0.001;
    out.injOXd = 0.001;
    if strcmp(scaling,'halved')
        out.F = 480/2; % DECISION
        out.injFUn = 1; % DECISION
        out.injOXn = 4; % DECISION
    elseif strcmp(scaling,'nominal')
        out.F = 480; % DECISION
        out.injFUn = 2; % DECISION
        out.injOXn = 8; % DECISION
    elseif strcmp(scaling,'doubled')
        out.F = 480*2; % DECISION
        out.injFUn = 4; % DECISION
        out.injOXn = 16; % DECISION
    end
end
% Gas constant [J/kg/K]
out.Rgas = R/out.Mm;
% Gamma recalculation, as the CEA value was said not to be reliable
gamma = 1/(1-out.Rgas/out.cp);

% Vandekerckhove and c*
out.Vandekerckhove = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)));
out.cstar = sqrt(out.Rgas*out.T(1))/out.Vandekerckhove;

% Solve isentropic flow from CC to nozzle Pa=0 boundary condition
options = optimset('Display','off');
f = @(y) epsPrelation(y, out.eps, gamma, out.Vandekerckhove);
PxPc = fsolve(f, 0.01, options);
out.P(4) = PxPc*out.P(1);
out.M(4) = sqrt(2/(gamma-1)*((1/PxPc)^((gamma-1)/gamma)-1));
out.M

% Find intermediate values for points A0, A, B, C
out.P(2:3) = out.P(1)./(1+(gamma-1)./2.*out.M(2:3).^2).^(gamma/(gamma-1));
out.T(2:4) = out.T(1).*(out.P(2:4)./out.P(1)).^((gamma-1)/gamma);
out.rho(2:4) = out.rho(1).*(out.P(2:4)./out.P(1)).^(1/gamma);
out.a = sqrt(gamma*out.Rgas.*out.T);
out.v = out.M.*out.a;

% Isp
out.Isp = (out.v(4)+out.P(4)/out.rho(4)/out.v(4))/g0;
% Mass flow for desired thrust and current Isp
out.mP = out.F/out.Isp/g0;

% Areas
out.A = [0;0;0;0];
% Throat area [m^2]
out.A(3) = out.mP/out.P(3)*sqrt(out.Rgas*out.T(3)/gamma);
% Exit area [m^2]
out.A(4) = out.A(3)*out.eps;
out.A(2) = out.A(3)*out.contraction;
out.A(1) = out.A(2);

% 2D effects
lambda2Dflow = (1+cosd(out.divergentAngle))/2;
out.F2D = lambda2Dflow*out.mP*out.v(4)+out.A(4)*out.P(4);
out.Isp2D = out.F2D/out.mP/g0;

out.gamma = gamma;
end

%% Geometry functions
function D = diamEngineConical(x, Lcc, Lcon, Dcc, Dt, convAngle, divAngle)
if x<=Lcc
    D = Dcc;
elseif x<=Lcc+Lcon
    D = Dcc-2*tand(convAngle)*(x-Lcc);
else
    D = Dt+2*tand(divAngle)*(x-Lcc-Lcon);
end
end

%% Integration functions
function res = solveEngine(A, x, Lmid, out)
if x<=Lmid
    y0 = 0.99;
else
    y0 = 0.01;
end

options = optimset('Display','off');
f = @(y) epsPrelation(y, A/out.A(3), out.gamma, out.Vandekerckhove);
PxPc = fsolve(f, y0, options);

P = PxPc*out.P(1);
M = sqrt(2/(out.gamma-1)*((1/PxPc)^((out.gamma-1)/out.gamma)-1));
T = out.T(1)*PxPc^((out.gamma-1)/out.gamma);
rho = out.rho(1)*PxPc^(1/out.gamma);
a = sqrt(out.gamma*out.Rgas*T);
v = M*a;

res = [M,P,T,rho,a,v];
end

function delta = epsPrelation(x, eps, gamma, Vandekerckhove)
delta = eps-Vandekerckhove/x^(1/gamma)/sqrt(2*gamma/(gamma-1)*(1-x^((gamma-1)/gamma)));
end

function delta = cWDarcyCoeff(f, D, Re, roughnessHeight)
fsqrt = sqrt(f);
delta = 1/fsqrt + 2*log10(roughnessHeight/3.7/D+2.51/Re/fsqrt);
end

%% Pressure dependent-heat capacities
function Tboil = TboilMMH(P)
% Clausius–Clapeyron equation: Tboil = f(T0, P0|T0, vapEnthalpy, P)
Tboil = 1/(1/364-8.314*log(P/100000)/(8.7621e+05));% [K]
end
function Tboil = TboilN2O4(P)
% Clausius–Clapeyron equation: Tboil = f(T0, P0|T0, vapEnthalpy, P)
Tboil = 1/(1/294.3-8.314*log(P/100000)/(4.3038e+05));% [K]
end
function Tboil = TboilRP1(P)
% Clausius–Clapeyron equation: Tboil = f(T0, P0|T0, vapEnthalpy, P)
Tboil = 1/(1/450-8.314*log(P/100000)/(2.2635e+05));% [K]
end
function Tcrit = TboilH2O2(P)
% This is a special topic. H2O2 decomposes into H2O and O2 at a certain
% rate. In a cool state, the decomposition results in liquid water and
% oxygen, but when we rise the T a boundary may be crossed, where the
% products of decomposition turn into boiling water, saturated steam
% and oxygen, which is undesirable. The following paper talks about it
% and gives a table of Tcritical vs Pressure.
% https://www.researchgate.net/publication/257615985_Thermodynamic_analysis_of_the_hydrogen_peroxide_decomposition_parameters/link/60e885451c28af3458595b53/download
Tcrit = 453.03+(485.52-453.03)*(P-10e6)/(20e6-10e6); % [K]
end

%% Temperature dependent-heat capacities
function cp = cpMMH(T)
% https://www.sciencedirect.com/science/article/pii/S0378381212001215?via%3Dihub
cp = 2832.6+(2930.3-2832.6)*(T-220)/(298-220); % [J/kg/K]
end
function cp = cpN2O4(T)
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C10544726&Mask=2&Type=JANAFL&Table=on
cp = 1548.7+(1943.2-1548.7)*(T-298)/(500-298); % [J/kg/K]
end
function cp = cpRP1(T)
% Heat capacity of rocket propellant (RP-1 fuel) at high temperatures and high pressures, I.M. Abdulagatov, N.D. Azizov, 2010
cp = 2016+(2298-2016)*(T-293.76)/(373.4-293.76); % [J/kg/K]
end
function cp = cpH2O2(T)
% Shomate coefficients https://booksite.elsevier.com/9780750683661/Appendix_C.pdf
cp = (-15.248+6.7693E-01*T-1.4948E-03*T^2+1.2018E-06*T^3)/(34.014/1000); % [J/kg/K]
end
function cp = cpGasMMH(T)
% https://hal.science/hal-03201924/document
cp = 2700;%1837.3+(2710.9-1837.3)*(T-500)/(1000-500); % [J/kg/K]
end
function cp = cpGasRP1(T)
% https://www.sciencedirect.com/science/article/pii/S0040603118301436
cp = 3500;%2200+(3500-2200)*(T-423)/(503-423); % [J/kg/K]
end

%% Plot functions
function plotOnEngine2D(dataStruc,titleStr,labelStr,xi,Di,D,Ltotal)
fn = fieldnames(dataStruc);

figure();
yyaxis left
if numel(fn)==1
    plot(xi/Ltotal,dataStruc.(fn{1}),'LineWidth',2,'LineStyle','-')
    hold on
else
    C1 = [0 0.4470 0.7410];  % blue
    C2 = [0.8500 0.3250 0.0980];  % orange
    C3 = [0.9290 0.6940 0.1250];  % yellow
    C4 = [0.4940 0.1840 0.5560];  % purple
    C5 = [0.4660 0.6740 0.1880];  % green
    C6 = [0.3010 0.7450 0.9330];  % light blue
    C7 = [0.6350 0.0780 0.1840];  % red
    for i=1:numel(fn)
        plot(xi/Ltotal,dataStruc.(fn{i}),'LineWidth',2,'DisplayName',fn{i},'LineStyle','-','Color',eval(['C' num2str(i)]))
        hold on
    end
end
title(titleStr,'Interpreter','latex');
xlabel('\boldmath$x/L$','Interpreter','latex'); ylabel(labelStr,'Interpreter','latex');
grid on;
yyaxis right
plot(xi/Ltotal,Di/D(4),'LineWidth',2,'Color','black','HandleVisibility', 'off')
hold on
plot([0,0],[0,D(1)/D(4)],'LineWidth',2,'Color','black','LineStyle','-','HandleVisibility', 'off')
hold on
plot([1,1],[0,1],'LineWidth',2,'Color','black','LineStyle','--','HandleVisibility', 'off')
hold on
plot([0,1],[0,0],'LineWidth',2,'Color','black','LineStyle','-.','HandleVisibility', 'off')
ylabel('\boldmath$y/R$','Interpreter','latex');
xlim([0 1])
if numel(fn)~=1
    legend;
end
set(gca,'fontsize', 16) 
hold off
end