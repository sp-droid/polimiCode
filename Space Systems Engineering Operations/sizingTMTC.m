%% Tracking Telemetry & Telecommand subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Sources
% S1: Detecting_Junos_Heartbeat_Communications_Support_during_Critical_Events_of_the_Juno_Mission.pdf
% S2: Juno Telecommunication.pdf
% S3: https://omniweb.gsfc.nasa.gov/coho/helios/heli.html
% S4: The_Juno_Gravity_Science_Instrument.pdf
% S5: https://en.wikipedia.org/wiki/Goldstone_Deep_Space_Communications_Complex#Antennas

%% Parameters
ANTENNA = "HGA";
LINK = "Uplink";
BAND = "X";
switch BAND
    case "X"
        switch LINK
            case "Uplink"
                FREQ_CARRIER = 7.153e9;             % Hz, S4
                DATARATE = 2e3;                     % bps, S2
            case "Downlink"
                FREQ_CARRIER = 8.404e9;
                DATARATE = 120e3;%18e3;
                % MOD INDEX SEEMS TO BE 0 IN Ka BAND!
        end
        AMPLIFIER = "TWTA";                         % S1
        AMP_OUT_POWER = 25;                         % W, S2
        AMP_MASS = 2.4;                             % kg, S2
        AMP_EFFICIENCY = 0.48;                      % Class notes
    case "Ka" % Unless specified, everything here has the same units / sources
        switch LINK
            case "Uplink"
                FREQ_CARRIER = 34.367e9;
                DATARATE = 2e3;
            case "Downlink"
                FREQ_CARRIER = 32.085e9;
                DATARATE = 18e3;
                % MOD INDEX SEEMS TO BE 0 IN Ka BAND!
        end
        AMPLIFIER = "SSPA";
        AMP_OUT_POWER = 2.5;
        AMP_MASS = 1.1;                             % Class notes
        AMP_EFFICIENCY = 0.11;
end

switch LINK
    case "Uplink"
        BER = 1e-5;                                 % Bit error rate, S2
        % This is the TRACKING LOOP BANDWIDTH, NOT THE NOISE BANDWIDTH!
        BANDWIDTH = 50;                             % Hz, S2
        CNR_MARGIN = 12;                            % dB
    case "Downlink"
        BER = 1e-6;
        BANDWIDTH = 50;
        CNR_MARGIN = 10;
end

MODULATION = "BPSK";                                % S2
MODULATION_COEF = 1;                                % Class notes
MODULATION_INDEX = 72;                              % º, S2
ENCODING = "Turbo 1/6, 8920 frame length";          % S2
% They do use other encoding schemes like concatenated convolutional + RS
% but it's for lower bit rates
ENCODING_COEF = 1/6;

switch ANTENNA % S2
    case "HGA"
        ANTENNA_D = 2.5;                            % deg, S2
        ANTENNA_BEAMWIDTH_REAL = 0.25;              % deg, S2
        ANTENNA_EFFICIENCY = 0.55; % Class notes for parabolic, even though this is a dual reflector
        ANTENNA_ACCURACY = 0.069;                   % deg, S2
        switch BAND
            case "X"
                ANTENNA_T = 400;                                 % K, S2
                switch LINK
                    case "Uplink"
                        ANTENNA_G_REAL = 43;
                    case "Downlink"
                        ANTENNA_G_REAL = 44.5;
                end
            case "Ka"
                ANTENNA_T = 700;
                switch LINK
                    case "Uplink"
                        ANTENNA_G_REAL = 47.5;
                    case "Downlink"
                        ANTENNA_G_REAL = 47;
                end
        end
    case "TLGA"
        ANTENNA_D = 0.26;
        ANTENNA_T = 400;
        ANTENNA_BEAMWIDTH_REAL = 10;
        ANTENNA_ACCURACY = 0.069;
        switch LINK
            case "Uplink"
                ANTENNA_G_REAL = 5.5;
            case "Downlink"
                ANTENNA_G_REAL = 6.5;
        end
    case "MGA"
        ANTENNA_D = 0.135;
        ANTENNA_EFFICIENCY = 0.52;
        ANTENNA_T = 400;
        ANTENNA_BEAMWIDTH_REAL = 10;
        ANTENNA_ACCURACY = 0.069;
        switch LINK
            case "Uplink"
                ANTENNA_G_REAL = 18.1;
                ANTENNA_BEAMWIDTH_REAL = 10.3;
            case "Downlink"
                ANTENNA_G_REAL = 18.8;
                ANTENNA_BEAMWIDTH_REAL = 9.3;
        end
    case "LGA"
        ANTENNA_D = 0.04;
        ANTENNA_EFFICIENCY = 0.52;
        ANTENNA_T = 400;
        ANTENNA_ACCURACY = 0.069;
        switch LINK
            case "Uplink"
                ANTENNA_G_REAL = 8.7;
                ANTENNA_BEAMWIDTH_REAL = 40;
            case "Downlink"
                ANTENNA_G_REAL = 7.7;
                ANTENNA_BEAMWIDTH_REAL = 42;
        end
end

DSN = "DSS Apollo 24, 25 & 26 (34m dish)";          % S5
DSN_D = 34;                                         % m, S5
DSN_EFFICIENCY = 0.55; % Class notes for parabolic
DSN_ACCURACY = 0.0055;

switch BAND
    case "X"
        DSN_OUT_POWER = 19999;                      % W, S2
        DSN_T = 33;                                 % K, S5
        switch LINK
            case "Uplink"
                DSN_G_REAL = 66.93;
            case "Downlink"
                DSN_G_REAL = 68.26;
        end
    case "Ka"
        DSN_OUT_POWER = 794;
        DSN_T = 95;
        switch LINK
            case "Uplink"
                DSN_G_REAL = 79.52;
            case "Downlink"
                DSN_G_REAL = 78.41;
        end
end

%% Loading the target coordinates S3
load('raw/earth.mat');
load('raw/jupiter.mat');
load('raw/juno.mat');

earth.x = earth.RAD_AU .* cosd(earth.HGI_LA) .* cosd(earth.THGI_LON);
earth.y = earth.RAD_AU .* cosd(earth.HGI_LA) .* sind(earth.THGI_LON);
earth.z = earth.RAD_AU .* sind(earth.HGI_LA);
jupiter.x = jupiter.RAD_AU .* cosd(jupiter.HGI_LA) .* cosd(jupiter.THGI_LON);
jupiter.y = jupiter.RAD_AU .* cosd(jupiter.HGI_LA) .* sind(jupiter.THGI_LON);
jupiter.z = jupiter.RAD_AU .* sind(jupiter.HGI_LA);
juno.x = juno.RAD_AU .* cosd(juno.HGI_LA) .* cosd(juno.THGI_LON);
juno.y = juno.RAD_AU .* cosd(juno.HGI_LA) .* sind(juno.THGI_LON);
juno.z = juno.RAD_AU .* sind(juno.HGI_LA);

time.year = juno.YEAR + juno.DAY/365;
distance.EarthToSun = earth.RAD_AU;
distance.JupiterToSun = jupiter.RAD_AU;
distance.JunoToSun = juno.RAD_AU;
distance.JunoToEarth = distanceObjects(earth, juno);
distance.JunoToJupiter = distanceObjects(jupiter, juno);

%% Constants
AU = 149597870700;%m
BOLTZMANN = 1.38e-23;%Ws/K

%% Sizing
fprintf(strcat("Link: ",LINK,", operating frequency: ",string(FREQ_CARRIER/1e9)," GHz, data rate: ",string(DATARATE/1e3)," kbps\n"))

% AMPLIFIER
input_power = AMP_OUT_POWER / AMP_EFFICIENCY;
fprintf(strcat("Amplifier (",AMPLIFIER,") input power / output power (efficiency): ",string(input_power)," / ", string(AMP_OUT_POWER)," W (",string(AMP_EFFICIENCY*100),"%%)\n"));
fprintf(strcat("Bit Error Rate: ",string(BER),"\n"))

% SYMBOL RATE
symbolRate = DATARATE * MODULATION_COEF / ENCODING_COEF;
fprintf(strcat("Symbol rate: ",string(symbolRate/1e3)," kSymbols per second\n"));

wavelength = 3e8/FREQ_CARRIER;

% JUNO ANTENNA

switch ANTENNA
    case "HGA"
        antennaGain = 10*log10(pi^2 * ANTENNA_D^2 * ANTENNA_EFFICIENCY / wavelength^2);
        antennaBeamwidth = 65.3 * wavelength / ANTENNA_D;
    case "TLGA"
        antennaR = ANTENNA_D / 2;
        antennaH = antennaR / 2;
        antennaGain = 5*log(antennaH/wavelength)+3.5;
        antennaBeamwidth = 72 * wavelength / ANTENNA_D;
    case "MGA"
        antennaGain = 10*log10(pi^2 * ANTENNA_D^2 * ANTENNA_EFFICIENCY / wavelength^2);
        antennaBeamwidth = 72 * wavelength / ANTENNA_D;
    case "LGA"
        antennaGain = 10*log10(pi^2 * ANTENNA_D^2 * ANTENNA_EFFICIENCY / wavelength^2);
        antennaBeamwidth = 72 * wavelength / ANTENNA_D;
end

fprintf(strcat("Antenna(",ANTENNA,") predicted / real gain (%%diff): ",string(antennaGain)," / ",string(ANTENNA_G_REAL)," dB (",string(antennaGain/ANTENNA_G_REAL*100-100),"%%)\n"));

fprintf(strcat("Antenna(",ANTENNA,") predicted / real beamwidth (%%diff): ",string(antennaBeamwidth)," / ",string(ANTENNA_BEAMWIDTH_REAL),"º (",string(antennaBeamwidth/ANTENNA_BEAMWIDTH_REAL*100-100),"%%)\n"));
% This big difference is due to us using the formulas of a simple parabolic antenna, when in reality the HGA is a dual reflector much more complex one

% GROUND STATION ANTENNA
dsnGain = 10*log10(pi^2 * DSN_D^2 * DSN_EFFICIENCY / wavelength^2);
fprintf(strcat("Ground station(",DSN,") predicted / real gain (%%diff): ",string(dsnGain)," / ",string(DSN_G_REAL)," dB (",string(dsnGain/DSN_G_REAL*100-100),"%%)\n"));
dsnBeamwidth = 65.3 * wavelength / DSN_D;
fprintf(strcat("Ground station(",DSN,") predicted beamwidth: ",string(dsnBeamwidth),"º\n"));

% POINTING LOSSES
% The pointing accuracy of the DSN antenna was calculated by taking the real pointing loss of -0.1 dB
% The one from Juno was fine tuned but mostly kept the same as the recommended 0.1º

lossesPointing = -12*(ANTENNA_ACCURACY/ANTENNA_BEAMWIDTH_REAL)^2;
lossesPointing = lossesPointing + -12*(DSN_ACCURACY/dsnBeamwidth)^2;
fprintf(strcat("Pointing losses: ",string(lossesPointing)," dB\n"))

% ATMOSPHERIC LOSSES (ignoring rain)
lossesAtm = -0.2; %This model is nice https://www.mathworks.com/help/phased/ref/gaspl.html, but we need to input height and temperature profiles...
fprintf(strcat("Atmospheric losses: ",string(lossesAtm)," dB\n"))

% CABLE LOSSES
lossesCables = -0.9;

% EIRP - Effective Isotropic Radiated Power & RECEIVER POWER
switch LINK
    case "Uplink"
        transmitterOutPower = DSN_OUT_POWER;
        transmitterGain = DSN_G_REAL;
        receiverGain = ANTENNA_G_REAL;
    case "Downlink"
        transmitterOutPower = AMP_OUT_POWER;
        transmitterGain = ANTENNA_G_REAL;
        receiverGain = DSN_G_REAL;
end
fprintf(strcat("Transmitter power: ",string(10*log10(transmitterOutPower)+30)," dBm\n"))
eirp = 10*log10(transmitterOutPower)+30 + transmitterGain + lossesCables;
fprintf(strcat("Effective Isotropic Radiated Power: ",string(eirp)," dBm\n"))

% FREE SPACE LOSSES - FROM THIS PART, THE DISTANCE MATTERS
fprintf(strcat("Max distance Juno-Earth: ",string(max(distance.JunoToEarth))," AU\n"))
maxDistance = max(distance.JunoToEarth)*AU;
lossesSpace = 20*log10(wavelength/4/pi/maxDistance);
fprintf(strcat("Space losses: ",string(lossesSpace)," dB\n"))

powerReceiver = eirp + receiverGain + lossesSpace + lossesPointing + lossesAtm;
fprintf(strcat("Receiver power: ",string(powerReceiver)," dBm\n"))

% SYSTEM NOISE DENSITY
switch LINK
    case "Uplink"
        N0receiver = 10*log10(BOLTZMANN*ANTENNA_T)+30;
    case "Downlink"
        N0receiver = 10*log10(BOLTZMANN*DSN_T)+30;
end
fprintf(strcat("Receiver system noise density: ",string(N0receiver)," dBm\n"));
fprintf(strcat("Pt/No: ",string(powerReceiver - N0receiver)," dB\n"));

% EBN0 - Error per Bit to Noise density
ebN0 = powerReceiver - N0receiver - 10*log10(DATARATE);
fprintf(strcat("Error per Bit to Noise density: ",string(ebN0)," dB\n"));
fprintf(strcat("Maximum bitrate for a minimum 3.39 dB (BER1e-6): ",string(10^((powerReceiver-3.39-N0receiver)/10)/1000),"kbps \n"));
fprintf(strcat("Maximum bitrate for a minimum 3.00 dB (BER1e-5): ",string(10^((powerReceiver-3-N0receiver)/10)/1000),"kbps \n"));

% CARRIER MODULATION INDEX REDUCTION
powerModLoss = 20*log10(cosd(MODULATION_INDEX));
if (BAND=="Ka")
    powerModLoss = 0;
end
fprintf(strcat("Carrier Modulation Index Reduction: ",string(powerModLoss)," dB\n"))

% CARRIER POWER
powerCarrier = powerReceiver + powerModLoss;
fprintf(strcat("Carrier power: ",string(powerCarrier)," dBm\n"))
% Pc/N0
powerCarrierNoiseRatio = powerCarrier-N0receiver;
fprintf(strcat("Pc/N0: ",string(powerCarrierNoiseRatio)," dB\n"))

% BANDWIDTH
fprintf(strcat("Bandwidth: ",string(BANDWIDTH)," Hz or ",string(10*log10(BANDWIDTH))," dB\n"))

% SIGNAL TO NOISE RATIO
snrCarrier = powerCarrier - N0receiver - 10*log10(BANDWIDTH);
fprintf(strcat("SNR: ",string(snrCarrier)," dB\n"))

snrMargin = snrCarrier - CNR_MARGIN;
fprintf(strcat("SNR margin > 3: ",string(snrMargin)," dB\n"))

%% Plots

performance.Distance = distance.JunoToEarth;
performance.FreeSpace = 20*log10(wavelength/4/pi./performance.Distance/AU);
performance.Pointing = ones(size(performance.Distance))*lossesPointing;
performance.Atmospheric = ones(size(performance.Distance))*lossesAtm;
performance.Circuit = ones(size(performance.Distance))*lossesCables;
% performance.Total = performance.FreeSpace + performance.Pointing + performance.Atmospheric + performance.Circuit;

performance2.Distance = distance.JunoToEarth;
Prx = eirp + receiverGain + performance.FreeSpace + performance.Pointing + performance.Atmospheric;
N0rx = N0receiver;
PtN0 = Prx - N0rx;
switch LINK
    case "Uplink"
        performance2.kbpsMaxforEbN0 = 10.^((PtN0-3)./10)/1000;
    case "Downlink"
        performance2.kbpsMaxforEbN0 = 10.^((PtN0-3.39)/10)/1000;
end
performance2.kbpsMaxforEbN0(performance2.kbpsMaxforEbN0 > 200) = 200;
performance3.Distance = distance.JunoToEarth;   
Pc = Prx + powerModLoss;
PcN0 = Pc - N0rx;
performance3.CNR = PcN0 - 10*log10(BANDWIDTH);
performance3.CNRthreshold = 3+ones(size(performance.Distance))*CNR_MARGIN;
cnrMargin = performance3.CNR - performance3.CNRthreshold;

plotOnEngine2D( ...
    time, ...
    performance, ...
    strcat("Losses for ",LINK,"-",string(FREQ_CARRIER*1e-9)," GHz (",BAND,"-band) ",ANTENNA," antenna",""), ...
    '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
    '\boldmath$Loss\hspace{0.5em}[dB]$', ...
    '\boldmath$Date\hspace{0.5em}[year]$');

% plotOnEngine2D( ...
%     time, ...
%     performance2, ...
%     strcat("Max bitrate for ",LINK,"-",string(FREQ_CARRIER*1e-9)," GHz (",BAND,"-band) ",ANTENNA," antenna",""), ...
%     '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
%     '\boldmath$Bitrate\hspace{0.5em}[kbps]$', ...
%     '\boldmath$Date\hspace{0.5em}[year]$');

% plotOnEngine2D( ...
%     time, ...
%     performance3, ...
%     strcat("Carrier loop SNR for ",LINK,"-",string(FREQ_CARRIER*1e-9)," GHz (",BAND,"-band) ",ANTENNA," antenna",""), ...
%     '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
%     '\boldmath$CNR\hspace{0.5em}[dB]$', ...
%     '\boldmath$Date\hspace{0.5em}[year]$');

%% Functions
function distance = distanceObjects(object1, object2)
    distance = sqrt((object1.x-object2.x).^2+(object1.y-object2.y).^2+(object1.z-object2.z).^2);
end

function plotOnEngine2D(graphx,graphy,titleStr,yLabelLeft,yLabelRight,xLabelStr)
fnx = fieldnames(graphx);
fny = fieldnames(graphy);

figure;
if numel(fny)==1
    scatter(graphx.(fnx{1}),graphy.(fny{1}))
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
    for i=1:numel(fny)
        if numel(fnx)==1
            fnxi = fnx{1};
        else
            fnxi = fnx{i};
        end
        if i==1
            yyaxis left
            plot(graphx.(fnxi),graphy.(fny{i}),'LineWidth',2,'DisplayName',fny{i},'Color',eval(['C' num2str(i)]))
            title(titleStr,'Interpreter','latex');
            xlabel(xLabelStr,'Interpreter','latex'); ylabel(yLabelLeft,'Interpreter','latex');
            yyaxis right
            continue
        end
        %scatter(graphx.(fnxi),graphy.(fny{i}),2,'DisplayName',fny{i},'MarkerEdgeColor',eval(['C' num2str(i)]))
        plot(graphx.(fnxi),graphy.(fny{i}),'LineWidth',2,'DisplayName',fny{i},'Color',eval(['C' num2str(i)]))
        hold on
    end
end
ylabel(yLabelRight,'Interpreter','latex');
grid on;
xlim([min(graphx.(fnx{1})), max(graphx.(fnx{1}))])

if numel(fny)~=1
    legend('Location','SouthEast');
end
set(gca,'fontsize', 16) 
hold off
end