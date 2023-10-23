clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Constants
% Real traces
load('tracesbar1.mat')
pExpTraces = who; % CAREFUL WITH THE WHO FUNCTION, IT RECORDS EVERY VAR IN THE WORKSPACE

% Fix trace mismatch
temp = pbar2444(:,3);
pbar2444(:,3) = pbar2444(:,2);
pbar2444(:,2) = temp;
clearvars temp

% Nominal pressure for traces
Pc = [30; 45; 70];

nTracesConstP = length(pExpTraces);
nPressures = min(size(eval(pExpTraces{1})));

% CEA and derivatives
R = 8.314; % [J/mol/K]
Mm = [25.574; 25.663; 25.758]/1000; % CEA
cp = [2.0456; 2.0486; 2.0519]*1000; % CEA
T = [3333.74; 3368.27; 3404.10];    % CEA

Rgas = R./Mm;
gamma = 1./(1-Rgas./cp);

Vandekerckhove = sqrt(gamma.*(2./(gamma+1)).^((gamma+1)./(gamma-1)));
cStarTheo = sqrt(Rgas.*T)./Vandekerckhove;

% Other parameters
Dext = 160; %mm
Dint = 100; %mm
L = 290; %mm
webThickness = 30; % mm
Dthroat = [28.8, 25.25, 21.81]; % mm
Athroat = pi/4*Dthroat.^2; % Circular
VburnTotal = L*pi/4*(Dext^2-Dint^2)/1e9; % Longitude 290mm, internal diam. 100mm, ext. diam. 160mm

density = 1/(0.68/1.95+0.18/2.7+0.14/0.92)*1000; % 68% AP, 18% Al, 14% HTPB binder
Mtotal = VburnTotal*density; % kg

%% Having loaded the experimental traces, go over each one of them and calculate P effective, tburn and Cstar as it is shown in the BC-V01.pdf 
tBurn = zeros(nTracesConstP,nPressures);
Peff = zeros(nTracesConstP,nPressures);
rb = zeros(nTracesConstP,nPressures);
cStarExp = zeros(nTracesConstP,nPressures);

pressureNames = {'LowP','MedP','HighP'};
for i = 1:length(pExpTraces)
    name = pExpTraces{i};
    traces = eval(name);

    time = (0:length(traces)-1)/1000';
    
    for j = 1:3
        tTraces.(pressureNames{j}).(pExpTraces{i}) = time;
        trace = traces(:,j);
        pTraces.(pressureNames{j}).(pExpTraces{i}) = trace-1.01325; % Relative pressure, minus 1 atm
    
        Pthreshold = 0.05*max(trace); % 5%
        idxA = 1+find(trace(1:round(length(trace)/2)) <= Pthreshold, 1, 'last');
        idxG = -1+round(length(trace)/2)+find(trace(round(length(trace)/2):end) <= Pthreshold, 1);

        tAction = time(idxG)-time(idxA);
        Pref = trapz(time(idxA:idxG), trace(idxA:idxG))/2/tAction;
        idxB = 1+find(trace(1:round(length(trace)/2)) <= Pref, 1, 'last');
        idxE = -1+round(length(trace)/2)+find(trace(round(length(trace)/2):end) <= Pref, 1);

        tBurn(i,j) = time(idxE)-time(idxB);
        Peff(i,j) = trapz(time(idxB:idxE), trace(idxB:idxE))/tBurn(i,j);
        rb(i,j) = webThickness/tBurn(i,j);

        cStarExp(i,j) = Peff(i,j)*1e5*tBurn(i,j)*Athroat(j)/1e6/Mtotal;
    end
end

%% Fit and uncertainty computation
cStar = mean(cStarExp(:));
etaCstar = cStar/cStarTheo(1);
Inc_cStar = std(cStarExp(:));
[a, Inc_a, n, Inc_n, R2] = Uncertainty(Peff(:), rb(:));
R2

% This is for a plot later
nPoints = 200;
Pplot = linspace(min(Peff(:)),max(Peff(:)),nPoints)';

rbPlot = a.*Pplot.^n;
rbUpper = (a+Inc_a).*Pplot.^(n+Inc_n);
rbLower = (a-Inc_a).*Pplot.^(n-Inc_n);

%% Monte-Carlo method
nPointsMonte = 30;
aMontes = normrnd(a, Inc_a^2, nPointsMonte, 1);
nMontes = normrnd(n, Inc_n^2, nPointsMonte, 1);
cStarMontes = normrnd(cStar, Inc_cStar^2, nPointsMonte, 1);

for j=1:3
    MEOPmc = zeros(nPointsMonte^3,1);
    MEOPmcMean = MEOPmc;
    MEOPmcInc = MEOPmc;
    deltat = 0.001;
    k=0;
    for iA=1:nPointsMonte
        aMonte = aMontes(iA);
        for iN=1:nPointsMonte
            nMonte = nMontes(iN);
            for iC=1:nPointsMonte
                k=k+1;
                cStarMonte = cStarMontes(iC);
                P = [];
                deltaBurnt = 0;
                i=0;
                while deltaBurnt<webThickness
                    i=i+1;
                    Ab = AbStep(Dext, Dint, L, deltaBurnt);
                    P(i) = (aMonte/1e8*density*Ab/Athroat(j)*cStarMonte)^(1/(1-nMonte));
                    rbi = aMonte*P(i)^(nMonte);
                    deltaBurnt = deltaBurnt + rbi*deltat;
                end
                %tBurnMC(k) = i*deltat;
                %tBurnMCiterMean(k) = mean(tBurnMC(1:k));
                %tBurnMCiterInc(k) = std(tBurnMC(1:k));
                MEOPmc(k) = max(P);
                MEOPmcMean(k) = mean(MEOPmc(1:k));
                MEOPmcInc(k) = std(MEOPmc(1:k));
            end
        end
    end
    
    if j==3
        figure;
        plot(1:k,MEOPmcMean,'Color','blue','LineWidth',1.5)
        title("\boldmath$"+num2str(Pc(j))+"\hspace{0.5em}bar\hspace{0.5em}\overline{P}_{max}="+num2str(MEOPmcMean(end))+"bar\hspace{0.5em}(Monte-Carlo\hspace{0.5em}method)$",'Interpreter','latex');
        xlabel('\boldmath$Iteration$','Interpreter','latex');
        ylabel('\boldmath$\overline{P}_{max} [bar]$','Interpreter','latex');
        grid on;
        set(gca,'fontsize', 16)
        
        figure;
        plot(1:k,MEOPmcInc,'Color','blue','LineWidth',1.5)
        title("\boldmath$"+num2str(Pc(j))+"\hspace{0.5em}bar\hspace{0.5em}\sigma_{P_{max}}="+num2str(MEOPmcInc(end))+"bar\hspace{0.5em}(Monte-Carlo\hspace{0.5em}method)$",'Interpreter','latex');
        xlabel('\boldmath$Iteration$','Interpreter','latex');
        ylabel('\boldmath$\sigma_{P_{max}} [bar]$','Interpreter','latex');
        grid on;
        set(gca,'fontsize', 16)
    end

    % Paint the nominal pressure traces with the now found burn times
    P = [0];
    time = [0];
    deltaBurnt = 0;
    i=1;
    while deltaBurnt<webThickness
        i=i+1;
        time(i) = time(i-1)+deltat;
        Ab = AbStep(Dext, Dint, L, deltaBurnt);
        %Vfree = (L-2*deltaBurnt)*2*pi*(Dint/2+deltaBurnt);
        P(i) = (a/1e8*density*Ab/Athroat(j)*cStar)^(1/(1-n));
        %P(i) = P(i-1) + deltat*(Ab*a*P(i-1)^n/Vfree*(density*Rgas(j)*T(j)-P(i-1)*1e5)-P(i-1)*1e5*Athroat(j)/Vfree*1e3*sqrt( ...
        %    gamma(j)*Rgas(j)*T(j)*(2/(gamma(j)+1))^((gamma(j)+1)/(gamma(j)-1))))*1e-5;
        rbi = a*P(i)^(n);
        deltaBurnt = deltaBurnt + rbi*deltat;
    end
    P(end+1) = 0; time(end+1) = time(end)+deltat;
    pTraces.(pressureNames{j}).Predicted = P;
    tTraces.(pressureNames{j}).Predicted = time;
end

figure;
scatter(Peff(:),rb(:),60,'DisplayName','Experimental data')
hold on
plot(Pplot, rbPlot, 'Color', 'black', 'LineStyle', '--','LineWidth',2,'DisplayName','Fit')
hold on
plot(Pplot, rbUpper, 'Color', 'black', 'LineStyle', '-.','LineWidth',0.1,'DisplayName','Fit upper bound')
hold on
fill([Pplot',fliplr(Pplot')],[rbPlot',fliplr(rbUpper')],'r','FaceAlpha',0.2,'HandleVisibility', 'off','EdgeAlpha',0)
hold on
plot(Pplot, rbLower, 'Color', 'black', 'LineStyle', '-.','LineWidth',0.1,'DisplayName','Fit lower bound')
hold on
fill([Pplot',fliplr(Pplot')],[rbPlot',fliplr(rbLower')],'r','FaceAlpha',0.2,'HandleVisibility', 'off','EdgeAlpha',0)
title('\boldmath$Pressure-burn\hspace{0.5em}rate$','Interpreter','latex');
xlabel('\boldmath$P_{eff}\hspace{0.5em}[bar]$','Interpreter','latex');
ylabel('\boldmath$r_{b}\hspace{0.5em}[\frac{mm}{s}]$','Interpreter','latex');
grid on;
legend('Location','SouthEast');
set(gca,'fontsize', 16)
hold off

plotOnEngine2D( ...
    tTraces.LowP, ...
    pTraces.LowP, ...
    '\boldmath$Time-Pressure\hspace{0.5em}trace\hspace{0.5em}(30\hspace{0.5em}bar)$', ...
    '\boldmath$P_{rel}\hspace{0.5em}[bar]$', ...
    '\boldmath$t\hspace{0.5em}[s]$');
plotOnEngine2D( ...
    tTraces.MedP, ...
    pTraces.MedP, ...
    '\boldmath$Time-Pressure\hspace{0.5em}trace\hspace{0.5em}(45\hspace{0.5em}bar)$', ...
    '\boldmath$P_{rel}\hspace{0.5em}[bar]$', ...
    '\boldmath$t\hspace{0.5em}[s]$');
plotOnEngine2D( ...
    tTraces.HighP, ...
    pTraces.HighP, ...
    '\boldmath$Time-Pressure\hspace{0.5em}trace\hspace{0.5em}(70\hspace{0.5em}bar)$', ...
    '\boldmath$P_{rel}\hspace{0.5em}[bar]$', ...
    '\boldmath$t\hspace{0.5em}[s]$');

%% Functions
function Ab = AbStep(Dext, Dint, L, deltaBurnt)
Ab = 2*pi*((Dext/2)^2-(Dint/2+deltaBurnt)^2) + (L-2*deltaBurnt)*2*pi*(Dint/2+deltaBurnt);
end

function plotOnEngine2D(graphx,graphy,titleStr,labelStr,xLabelStr)
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
        if i==10
            plot(graphx.(fn{i}),graphy.(fn{i}),'LineWidth',2,'DisplayName',fn{i},'LineStyle','-','Color',eval(['C' num2str(i)]))
        else
            scatter(graphx.(fn{i}),graphy.(fn{i}),2,'DisplayName',fn{i},'MarkerEdgeColor',eval(['C' num2str(i)]))
        end
        hold on
    end
end
title(titleStr,'Interpreter','latex');
xlabel(xLabelStr,'Interpreter','latex'); ylabel(labelStr,'Interpreter','latex');
grid on;

if numel(fn)~=1
    legend('Location','South');
end
set(gca,'fontsize', 16) 
hold off
end