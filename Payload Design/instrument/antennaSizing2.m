%% Antenna sizing
clear
close all
clc
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

%% Constants
c = 2.998E+8;                               % Speed of light in vacuum
muEarth = astroConstants(13);               % Earth gravity constant
Rearth = astroConstants(23)*1E3;            % Earth radius
k = 1.380649E-23;                           % Boltzmann constant

%% Parameters
thetaFar = 3.915;                           % [deg] Far range look angle
H = 892.09e3;                               % [m] Orbital height
vOrb = sqrt(muEarth/(Rearth+H)/1e-3)*1e3;   % [m/s] Orbital velocity

nPolarizations = 2;                         % [-] No. polarizations. VV and HH for the two sides
thetaNear = 0.65;                           % [deg] Near range look angle
ReAT = 2.5;                                 % [m/pixel] Along-track resolution
ReCT = 60;                                  % [m/pixel] Cross-track resolution
targetSNR = 17;                             % [dB] Target signal-to-noise ratio

effAperture = 0.5;                          % [-] Aperture efficiency
TxB = 1100;                                 % [-] Chirp compression factor
Tsys = 290;                                 % [K] Scene temperature
ka = 1.5;                                   % [-] PRF min. aliasing correction factor
F = 6;                                      % [dB] System noise factor
Lgrange = 1.3;                              % [dB] Non-ideal range filtering loss
Lrain = 5;                                  % [dB] Rain attenuation
Lradar = 2;
scatterDB = 15;

pulseFreq = 33E9;

%%
L = 2*ReAT;
fprintf(strcat("Antenna length: ",string(L)," m\n"))

B = c/2/ReCT/sind(thetaNear);
fprintf(strcat("Bandwidth for ",string(ReCT)," m near-end resolution: ",string(B*1e-6)," MHz\n"))
Ntiles = round(2*B*H/c*(cosd(thetaNear)-cosd(thetaFar)));
fprintf(strcat("Cross-track resolution (near and far end): ",string(ReCT)," -- ",string(c/2/B/sind(thetaFar))," m, in ",string(Ntiles)," tiles\n"))
tPulse = TxB/B;

swathReachFar = H*tand(thetaFar);
swathReachNear = H*tand(thetaNear);
SW = swathReachFar-swathReachNear;
fprintf(strcat("Swath width: ",string(swathReachNear*1e-3)," <--",string(SW*1e-3),"--> ",string(swathReachFar*1e-3)," km\n"))

pulseWavelength = c/pulseFreq;

W = pulseWavelength/deg2rad(thetaFar-thetaNear)/effAperture;
fprintf(strcat("Antenna width: ",string(W)," m\n"))

antennaMass = 20.5044*L*W;
fprintf(strcat("Antenna mass: ",string(antennaMass)," kg\n"))

PRFmin = ka*vOrb/ReAT;
PRFmax = 1/(tPulse+2*H/c*(1/cosd(thetaFar)-1/cosd(thetaNear)));
PRFmaxUnambiguous = 1/(tPulse+2*H/c/cosd(thetaFar));
fprintf(strcat("PRF min: ",string(PRFmin)," Hz, max: ",string(PRFmax)," Hz, max unambiguous: ",string(PRFmaxUnambiguous)," Hz\n"))
PRF = PRFmin;
duty = PRF*tPulse;
fprintf(strcat("Duty cycle: ",string(duty*100)," %%\n"))

range = H/cosd(thetaFar);

NintegratedPulses = PRF*pulseWavelength*range/2/ReAT/vOrb;
tIntegration = NintegratedPulses/PRF;
Gazimuth = NintegratedPulses;
fprintf(strcat("Number of integrated pulses (G_a): ",string(NintegratedPulses)," (",string(10*log10(NintegratedPulses))," dB), Integration time: ",string(tIntegration)," s\n"))

Lgrange = 10^(Lgrange/10);
Grange = TxB/Lgrange;
fprintf(strcat("Range processing gain: ",string(Grange)," (",string(10*log10(Grange))," dB)\n"))

Q = (8+16);
datarate = 2*2*B*32*1E-9*2*PRF*duty*2*SW/c;%x2 antennas receiving, emission at x2 PRF  %Ntiles*Q*2*2*PRF*1E-9;
datarateLR = datarate/310*1000;
datarateHR = datarate/13*1000;
fprintf(strcat("Datarate: ",string(datarate)," gbps, LR datarate: ",string(datarateLR)," mbps, max unambiguous: ",string(datarateHR)," mbps\n"))

F = 10^(F/10);
Pnoise = k*F*Tsys*B;
Psignal = 10^(targetSNR/20)*Pnoise/Gazimuth/Grange;

G = effAperture*4*pi*L*W/pulseWavelength^2;
scatter = 10^(scatterDB/10);

Latm = 2*10^(pathAttenuationITUR676(90 - thetaFar, pulseFreq*1E-9, 15, 101325, 20)/10+Lrain/10);
Lradar = 10^(Lradar/10);

Ppeak = Psignal*(4*pi)^3*range^4*Latm*Lradar/scatter/G^2/pulseWavelength^2;
Pavg = Ppeak*nPolarizations*duty;

fprintf(strcat("Peak transmit power: ",string(Ppeak)," W, Average transmit power: ",string(Pavg)," W\n"))

% % Baseline vs height error
% baselines = linspace(1,12,200);
% deltaHs = zeros(1,length(baselines));
% for i=1:length(baselines)
%     baseline = baselines(i);
%     theta = thetaFar;
%     C = H*tand(theta);
% 
%     SWH = 2;
%     SWHstd = SWH/4;
%     phiSTDwave = 2*pi*SWHstd*baseline*cosd(theta)/H/tand(theta)/pulseWavelength;
%     gammaWave = exp(-0.5*phiSTDwave^2);
%     SNR = 10^(17/20);
%     gammaThermal = 1/(1+1/SNR);
%     gammaTime = 0.8;
%     gammaBaseline = 1-pulseFreq*baseline*cosd(theta)^2/(B*H*tand(theta));
%     gamma = gammaWave*gammaThermal*gammaTime*gammaBaseline;
%     phiSTD = sqrt((1-gamma^2)/2/gamma^2/24/Ntiles);
%     deltaHs(i) = phiSTD*C/2/pi*pulseWavelength/baseline*(1+H/Rearth);
% end
% 
% figure;
% plot(baselines, deltaHs*100, 'LineWidth', 2)
% hold on;
% plot(baselines, 3*ones(1,length(baselines)))
% grid on;
% xlabel('Baseline length [m]')
% ylabel('Height resolution [cm]')
% xlim([min(baselines),max(baselines)])
% set(gca, 'fontsize', 14)
% print(gcf, 'output/baselineHeightError.png', '-dpng', '-r200');
% 
% Height error vs swath
thetas = linspace(thetaNear,thetaFar,200);
swath = H*tand(thetas)*1E-3;
deltaHs = zeros(1,length(thetas));
for i=1:length(baselines)
    baseline = 10;
    theta = thetas(i);
    C = H*tand(theta);

    SWH = 2;
    SWHstd = SWH/4;
    phiSTDwave = 2*pi*SWHstd*baseline*cosd(theta)/H/tand(theta)/pulseWavelength;
    gammaWave = exp(-0.5*phiSTDwave^2);
    SNR = 10^(17/20);
    gammaThermal = 1/(1+1/SNR);
    gammaTime = 0.8;
    gammaBaseline = 1-pulseFreq*baseline*cosd(theta)^2/(B*H*tand(theta));
    gamma = gammaWave*gammaThermal*gammaTime*gammaBaseline;
    phiSTD = sqrt((1-gamma^2)/2/gamma^2/24/Ntiles);
    deltaHs(i) = phiSTD*C/2/pi*pulseWavelength/baseline*(1+H/Rearth);
end

figure;
plot(swath, deltaHs*100, 'LineWidth', 2)
hold on;
% plot(swath, 3*ones(1,length(swath)))
grid on;
xlabel('Swath [km]')
ylabel('Height resolution [cm]')
xlim([9,61])
% set(gca, 'fontsize', 14)
% print(gcf, 'output/SwathHeightError.png', '-dpng', '-r200');

%%
% freqs = linspace(20,50,200);
% powers = zeros(1,length(freqs));
% for i=1:length(freqs)
%     pulseWavelength = c/freqs(i)/1E9;
% 
%     W = pulseWavelength/deg2rad(thetaFar-thetaNear)/effAperture;
% 
%     antennaMass = 20.5044*L*W;
% 
%     PRFmin = ka*vOrb/ReAT;
%     PRFmax = 1/(tPulse+2*H/c*(1/cosd(thetaFar)-1/cosd(thetaNear)));
%     PRFmaxUnambiguous = 1/(tPulse+2*H/c/cosd(thetaFar));
%     PRF = PRFmin;
%     duty = PRF*tPulse;
% 
%     range = H/cosd(thetaFar);
% 
%     NintegratedPulses = PRF*pulseWavelength*range/2/ReAT/vOrb;
%     tIntegration = NintegratedPulses/PRF;
%     Gazimuth = NintegratedPulses;
% 
%     Grange = TxB/Lgrange;
% 
%     Q = (8+8);
%     datarate = Ntiles*Q*2*2*PRF*1E-9;
%     datarateLR = datarate/310*1000;
%     datarateHR = datarate/13*1000;
% 
%     Pnoise = k*F*Tsys*B;
%     Psignal = 10^(targetSNR/20)*Pnoise/Gazimuth/Grange;
% 
%     G = effAperture*4*pi*L*W/pulseWavelength^2;
% 
%     Latm = 2*10^(pathAttenuationITUR676(90 - thetaFar, freqs(i), 15, 101325, 20)/10+Lrain/10);
% 
%     Ppeak = Psignal*(4*pi)^3*range^4*Latm*Lradar/scatter/G^2/pulseWavelength^2;
%     powers(i) = Ppeak*nPolarizations*duty;
% 
% end
% 
% figure;
% plot(freqs, powers, 'LineWidth', 2)
% grid on;
% xlabel('Frequency [GHz]')
% ylabel('Average transmitted power [W]')
% set(gca, 'fontsize', 14)
% print(gcf, 'output/freqvsPower.png', '-dpng', '-r200');

%% execute this only after commenting the previous section IMPORTANT
thetas = linspace(thetaNear,thetaFar,200);
swath = H*tand(thetas)*1E-3;
SNRs = zeros(1,length(thetas));
rangeRes = zeros(1,length(thetas));
for i=1:length(thetas)
    theta = thetas(i);

    pulseWavelength = c/pulseFreq;

    range = H/cosd(theta);

    NintegratedPulses = PRF*pulseWavelength*range/2/ReAT/vOrb;
    tIntegration = NintegratedPulses/PRF;
    Gazimuth = NintegratedPulses;

    Grange = TxB/Lgrange;

    G = effAperture*4*pi*L*W/pulseWavelength^2;

    Latm = 2*10^(pathAttenuationITUR676(90 - theta, pulseFreq*1E-9, 15, 101325, 20)/10+Lrain/10);

    Psignal = Ppeak*1/((4*pi)^3*range^4*Latm*Lradar/scatter/G^2/pulseWavelength^2);
    Pnoise = k*F*Tsys*B;
    SNRs(i) = 20*log10(Psignal/Pnoise*Gazimuth*Grange);

    rangeRes(i) = c/2/B/sind(theta);
end

figure;
plot(swath, SNRs, 'LineWidth', 2)
grid on;
xlabel('Swath [km]')
ylabel('SNR [dB]')
xlim([9,61])
% set(gca, 'fontsize', 14)
% print(gcf, 'output/swathSNR.png', '-dpng', '-r200');
% 
% figure;
% plot(swath, rangeRes, 'LineWidth', 2)
% grid on;
% xlabel('Swath [km]')
% ylabel('Range resolution [m]')
% xlim([9,61])
% set(gca, 'fontsize', 14)
% print(gcf, 'output/swathResolution.png', '-dpng', '-r200');