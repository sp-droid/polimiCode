clear;
clc;
close all;
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% ENVIRONMENTAL DATA
g0 = 9.80665;
Tearth = 23*3600+56*60+4;
Rearth = 6378;

%% INPUTS
Isp = {[263,302], [199,419], [202,426]};              % Specific impulse
eps = {0.059, 0.08, 0.111};         % Inert mass fraction
mPL = 50000;                           % Payload mass in tons
vOrb = 8.728;                       % Target ideal delta V in km/s

%launchLat = 28;                     % Launch site latitude in º (Canary Islands)
h = {0,40,100};

%%
lambdas = linspace(0,2,10000);
deltaVs = arrayfun(@(x) vOrb-getdeltaV(x, Isp, eps, g0, mPL, h), lambdas);
[minValue, minIndex] = min(abs(deltaVs));
minValue
lambda = lambdas(minIndex)
%%
disp('--------')

c = cellfun(@(i) getc_i(h{i}, Isp{i}, g0), num2cell(1:length(eps)), 'UniformOutput', false);
n = cellfun(@(i) getn_i(lambda, c{i}, eps{i}), num2cell(1:length(eps)), 'UniformOutput', false);
m = stagesMass(n, eps, mPL);
m0 = sum(cellfun(@double, m))+mPL;

deltaV1 = c{1}*log(n{1});
[h1,h1Mean,~] = burnout1(m0, m{1}, eps{1}, c{1}, g0);
[h2,h2Mean,~] = burnout2(h1,m0-m{1}, m{2}, eps{2}, c{2}, g0, deltaV1);

h1;
h2;
hh = {h1Mean, h2Mean, h2};

c = cellfun(@(i) getc_i(hh{i}, Isp{i}, g0), num2cell(1:length(eps)), 'UniformOutput', false);
n = cellfun(@(i) getn_i(lambda, c{i}, eps{i}), num2cell(1:length(eps)), 'UniformOutput', false);
m = stagesMass(n, eps, mPL);
m0 = sum(cellfun(@double, m))+mPL;

deltaV = cellfun(@(i) c{i}*log(n{i}), num2cell(1:length(eps)), 'UniformOutput', false);
deltaVcheck = sum(cellfun(@double, deltaV));

%% ROBUST DESIGN METHODOLOGY
deltaVloss = 0.93;
deltaVideal = vOrb-deltaVloss;

[alpha, beta] = meshgrid(0:0.005:1, 0:0.005:1);
mOptim = alpha;

staging1 = [];
staging1m = [];
staging2 = [];
staging2m = [];
k1 = 0;
k2 = 0;

for i=1:length(alpha)
    for j=1:length(beta)
        deltaV = {deltaVideal*alpha(i,j)+deltaVloss*beta(i,j), deltaVideal*(1-alpha(i,j))+deltaVloss*(1-beta(i,j))};
        c = cellfun(@(i) getc_i(h{i}, Isp{i}, g0), num2cell(1:2), 'UniformOutput', false);
        n = cellfun(@(i) exp(deltaV{i}/c{i}), num2cell(1:2), 'UniformOutput', false);
        m = stagesMass(n, eps, mPL);
        m0 = sum(cellfun(@double, m))+mPL;
        
        [h1,h1Mean,~] = burnout1(m0, m{1}, eps{1}, c{1}, g0);
        [h2,h2Mean,~] = burnout2(h1,m0-m{1}, m{2}, eps{2}, c{2}, g0, deltaV{1});
        hh = {h1Mean, h2Mean, h2};
        
        c = cellfun(@(i) getc_i(hh{i}, Isp{i}, g0), num2cell(1:2), 'UniformOutput', false);
        n = cellfun(@(i) exp(deltaV{i}/c{i}), num2cell(1:2), 'UniformOutput', false);
        m = stagesMass(n, eps, mPL);
        if min(cellfun(@double,m))<0
            mOptim(i,j) = nan;
            continue
        end
        m0 = sum(cellfun(@double, m))+mPL;
        mOptim(i,j) = m0/1000;
        if mOptim(i,j) > 1500
            mOptim(i,j) = nan;
        end

        if beta(i,j)==0
            k1 = k1 + 1;
            staging1m(k1) = mOptim(i,j);
            staging1(k1) = deltaVideal*alpha(i,j);
        end
        if beta(i,j)==alpha(i,j)
            k2 = k2 + 1;
            staging2m(k2) = mOptim(i,j);
            staging2(k2) = deltaVideal*alpha(i,j);
        end
    end
end
[minValue, linearIndex] = min(mOptim(:));
[rowIndex, colIndex] = ind2sub(size(mOptim), linearIndex);
alphaOptim = alpha(rowIndex, colIndex)
betaOptim = beta(rowIndex, colIndex)
m0Optim = mOptim(rowIndex, colIndex)

figure;
contourf(alpha, beta, mOptim);
hold on;
plot(alphaOptim, betaOptim, 'ro', 'MarkerSize', 20);

xlabel('Alpha')
ylabel('Beta')

%%
figure;
plot(staging1,staging1m,'LineWidth',2,'DisplayName','beta = 0')
hold on;
plot(staging2,staging2m,':','LineWidth',2,'DisplayName','alpha = beta')
xlabel('Staging speed [km/s]')
ylabel('Take-off mass [ton]')
set(gca,'fontsize', 15)
legend;
grid on

%%
figure;
plot(lambdas,deltaVs,'LineWidth',2)
xlabel('Lambda')
ylabel('DeltaV-target')
set(gca,'fontsize', 15)
grid on
%%
wEarth = 2*pi/Tearth;
%cosd(launchLat);
%vEquator = wEarth*Rearth*cosd(launchLat);

%% FUNCTIONS
function Isp = getIsp_i(h, Isp_sl, Isp_vac)
    if h<100
        Isp = Isp_sl*exp(h/100*log(Isp_vac/Isp_sl));
    else
        Isp = Isp_vac;
    end
end

function c_i = getc_i(h, Isp_i, g0)
    Isp = getIsp_i(h, Isp_i(1), Isp_i(2));
    c_i = Isp*g0/1000;
end

function m_i = getm_i(n_i, eps_i, mUpper)
    m_i = (n_i-1)/(1-n_i*eps_i)*mUpper;
end

function m = stagesMass(n, eps, mPL)
    m = {};
    for i=length(n):-1:1
        mUpper = mPL;
        for j=length(n):-1:i+1
            mUpper = mUpper + m{j};
        end
        m{i} = getm_i(n{i},eps{i},mUpper);
    end
end

function n_i = getn_i(lambda, c_i, eps_i)
    n_i = (c_i*lambda-1)/c_i/eps_i/lambda;
end

function [h,hMean,acc] = burnout1(m0, m1, eps1, c1, g0)
    T1 = 34500000;
    mBurnt = m1*(1-eps1);
    t1 = mBurnt*c1*1000/T1;
    TWR = T1/m0/g0;
    acc = (TWR-1)*g0;
    h = 0.5*acc*t1^2/1000; % Height after first stage

    hMean = integral(@(t) 0.5*acc*t.^2/1000, 0, t1)/t1;
end

function [h,hMean,acc] = burnout2(h1, m01, m2, eps2, c2, g0, deltaV1)
    T2 = 4900000;
    mBurnt = m2*(1-eps2);
    t2 = mBurnt*c2*1000/T2;
    TWR = T2/m01/g0;
    acc = (TWR-1)*g0;
    h = h1 + deltaV1*t2 + 0.5*acc*t2^2/1000; % Height after second stage

    hMean = integral(@(t) h1 + deltaV1*t + 0.5*acc*t.^2/1000, 0, t2)/t2;
end

function deltaV = getdeltaV(lambda, Isp, eps, g0, mPL, h)
    deltaV = nan;
    if lambda>10 || lambda<0
        return
    end
    %h = {0,100,100};

    c = cellfun(@(i) getc_i(h{i}, Isp{i}, g0), num2cell(1:length(eps)), 'UniformOutput', false);
    if lambda < 1/min(cellfun(@double, c))
        return
    end
    n = cellfun(@(i) getn_i(lambda, c{i}, eps{i}), num2cell(1:length(eps)), 'UniformOutput', false);
    m = stagesMass(n, eps, mPL);
    if min(cellfun(@double, m))<0
        return % Negative masses
    end
    m0 = sum(cellfun(@double, m))+mPL;
    
    deltaV1 = c{1}*log(n{1});
    [h1,h1Mean,acc] = burnout1(m0, m{1}, eps{1}, c{1}, g0);
    if acc<0
        return % Thrust to weight ratio < 1
    end
    [h2,h2Mean,acc] = burnout2(h1,m0-m{1}, m{2}, eps{2}, c{2}, g0, deltaV1);
    if acc<0
        return % Thrust to weight ratio < 1
    end
    h{1} = h1Mean;
    h{2} = h2Mean;
    h{3} = h2;
    
    deltaV = 0;
    for i=1:length(eps)
        c = getc_i(h{i}, Isp{i}, g0);
        n = getn_i(lambda, c, eps{i});
        if lambda < 1/c
            deltaV = nan;
            return
        end
        deltaV = deltaV + c*log(n);
    end
    
end