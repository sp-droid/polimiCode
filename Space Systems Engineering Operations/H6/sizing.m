clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Load New Horizons RTG output dataset
load realRTG.dat

%% Mission constraints
powerBOM = 245.7;       % W
timeEOM = 9.5;          % Years
powerEOM = 201.7;         % W
powerEOMminimum = 191;  % W
time2022 = 16;          % years

%% Power source selection:
% RTG is selected as energy source over solar arrays because of the extreme
% distance to the Sun. RTG's output is not intermitent -> No need for
% batteries unless we need very high power peaks (we do not)

%% RTG characteristics:
% Source: class notes
% GE-RTG (General Electric RTG): m = 56 kg, P = 4400 W, eta = ~0.068
% MM-RTG (Multi Mission RTG): m = 43 kg, P = 2000 W, eta = ~0.063

% Both of them use 7.8 kg of Plutonium (Pu-238) as radioactive source, which as 87.7 of half-life decay period.
tau = 87.74;                  % Pu-238 half-life in years
lambda = log(2)/tau;         % decay constant

% Source: http://large.stanford.edu/courses/2015/ph241/chirayath1/docs/aiaa-2006-4029.pdf
% "F-8" GPHS-RTG
% 10.9 kg Pu-238 dioxide, ended up being 9.75 kg due to fuel processing
% delays. 15% of initial baseline power lost
% designed and built by the Johns Hopkins University Applied Physics Laboratory (APL)
% Annual loss of 0.8%
annualLoss = (1-exp(-lambda*1))*100 % This is 0.79%, which complies with the reference of 0.8%

eta = 0.064;                % Conversion efficiency
Q0 = 3948;                  % W
mass = 10.9;%9.75;                % kg
%mass = 7.8; % From the curie number 132465
% Each module contains 4 pellets, in total about 0.54 kg of Pu-238 dioxide
nModules = 18;
moduleMass = 0.605;%0.542;         % kg. Obtained from Q0/nModules/powerDensity
powerDensity = 405;         % W/kg

%% Visualization
nPoints = 200;
time = linspace(0,timeEOM,nPoints);

% Reality
Ptraces.Real = realRTG(:,2)*powerBOM;
tTraces.Real = realRTG(:,1);
disp('PowerBOM and PowerEOM, reality')
disp(realRTG(1,2)*powerBOM)
disp(realRTG(end,2)*powerBOM)

% Simplified model
P = 0.06223*Q0*exp(-lambda*time);
Ptraces.SimpleModel = P;
tTraces.SimpleModel = time;
disp('PowerBOM and PowerEOM, simple model')
disp(P(1))
disp(P(end))

% Advanced model, accounting for thermoelectric deg. and better fuel deg.
% https://www.osti.gov/servlets/purl/1033353
% This one-> https://www.osti.gov/servlets/purl/1047828
% Natural degradation and power conversion: Pu = aP + bP*Q0*exp(-lambda*t)
aP = -222.3;                % W
%bP = 0.1192;
bP = 0.1186; % The original coef. is for fully unobstructed, this may be better
% Instantaneous thermoelectric degradation: ratio = 1 - alpha*sqrt(t)
% Proportionality factor alpha = aAlpha*exp(-bAlpha/T)
aAlpha = 6973;              % years^0.5
bAlpha = 15480;             % kelvin
% Average hot-junction temperature T = aT + bT*Q0*exp(-lambda*t)
aT = 566.25;                % kelvin
bT = 0.1572;                % kelvin/W
% Thermoelectric degradation: ratio = 1 - sqrt(integral(alpha^2))
fun = @(t) exp(-2*bAlpha./(aT+bT*Q0*exp(-lambda*t)));
fun = @(t) integral(fun,0,t);

P = (aP+bP.*Q0.*exp(-lambda.*time)).*(1-aAlpha.*sqrt(arrayfun(fun, time)));
Ptraces.MoweryModel = P;
tTraces.MoweryModel = time;
disp('PowerBOM and PowerEOM, 2nd model')
disp(P(1))
disp(P(end))

diffBOM = (P(1)-powerBOM)/powerBOM*100
diffEOM = (P(end)-powerEOM)/powerEOM*100

P = (aP+bP.*Q0.*exp(-lambda.*time2022)).*(1-aAlpha.*sqrt(fun(time2022)));
disp('Power2022, 2nd model')
disp(P)

plotSpecial( ...
    tTraces, ...
    Ptraces, ...
    "\boldmath$RTG's\hspace{0.5em}power\hspace{0.5em}degradation$", ...
    '\boldmath$Electrical\hspace{0.5em}power\hspace{0.5em}output\hspace{0.5em}[W]$', ...
    '\boldmath$Time\hspace{0.5em}after\hspace{0.5em}launch\hspace{0.5em}[years]$');

%% Module sizing using the EOM power requirement
fun = @(estQ0,t) exp(-2*bAlpha./(aT+bT*estQ0*exp(-lambda*t)));
fun = @(estQ0) integral(@(t) fun(estQ0,t),0,timeEOM);
fun = @(estQ0) -(201)+(aP+bP.*estQ0.*exp(-lambda.*timeEOM)).*(1-aAlpha.*sqrt(fun(estQ0)));
disp('---')
estQ0 = fsolve(fun, 4000)
diffQ0 = (estQ0-Q0)/Q0*100

% Mass considering power density
disp('---')
estMass = estQ0/powerDensity
diffMass = (estMass-mass)/mass*100

disp('---')
estNmodules = ceil(estMass/moduleMass)
diffModules = nModules-estNmodules

% Mass if all modules have specific mass and are considered equal
disp('---')
estMass = estNmodules*moduleMass
diffMass = (estMass-mass)/mass*100

%% Functions
function plotSpecial(graphx,graphy,titleStr,labelStr,xLabelStr)
fn = fieldnames(graphx);

figure;
if numel(fn)==1
    %plot(xi,graphy.(fn{1}),'LineWidth',2,'LineStyle','-')
    scatter(graphx.(fn{1}),graphy.(fn{1}))
    hold on
else
    C1 = [0 0.4470 0.7410];  % blue
    C2 = [0.8500 0.3250 0.0980];  % orange
    C3 = [0.9290 0.6940 0.1250];  % yellow
    C4 = [0.4940 0.1840 0.5560];  % purple
    C5 = [0.4660 0.6740 0.1880];  % green
    C6 = [0.3010 0.7450 0.9330];  % light blue
    C7 = [0.6350 0.0780 0.1840];  % red
    C8 = [1.0000 0.4000 0.6000];  % pink
    C9 = [0.6350 0.0780 0.1840];  % brown
    C10 = [0.0000 0.0000 0.0000];  % black
    for i=1:numel(fn)
        if i~=1
            plot(graphx.(fn{i}),graphy.(fn{i}),'LineWidth',2,'DisplayName',fn{i},'LineStyle','-','Color',eval(['C' num2str(i)]))
        else
            scatter(graphx.(fn{i}),graphy.(fn{i}),2,'filled','DisplayName',fn{i},'MarkerEdgeColor',eval(['C' num2str(i)]))
        end
        hold on
    end
end
title(titleStr,'Interpreter','latex');
xlabel(xLabelStr,'Interpreter','latex'); ylabel(labelStr,'Interpreter','latex');
grid on;

if numel(fn)~=1
    legend('Location','NorthEast');
end
set(gca,'fontsize', 16) 
hold off
end



