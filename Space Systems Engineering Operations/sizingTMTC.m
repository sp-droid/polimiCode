%% Tracking Telemetry & Telecommand subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.25 0.25 0.3 0.4]);
addpath('..\..\commonFunctions')

%% Sources
% S1: Detecting_Junos_Heartbeat_Communications_Support_during_Critical_Events_of_the_Juno_Mission.pdf
% S2: Juno Telecommunication.pdf
% S3: https://omniweb.gsfc.nasa.gov/coho/helios/heli.html
% S4: The_Juno_Gravity_Science_Instrument.pdf
% S5: https://en.wikipedia.org/wiki/Goldstone_Deep_Space_Communications_Complex#Antennas

%% Parameters
ANTENNA = "HGA";
LINK = "Downlink";
BAND = "X";
switch BAND
    case "X"
        switch LINK
            case "Uplink"
                FREQ_CARRIER = 7.153e9;             % Hz, S4
                DATARATE = 2e3;                     % bps, S2
            case "Downlink"
                FREQ_CARRIER = 8.404e9;
                DATARATE = 18e3;%18e3;
                % MOD INDEX SEEMS TO BE 0 IN Ka BAND!
        end
        AMPLIFIER = "TWTA";                         % S1
        AMP_OUT_POWER = 25;                         % W, S2
        AMP_MASS = 2.4;                             % kg, S2
        AMP_EFFICIENCY = 0.48;                      % Class notes
end

BER = 1e-5;                                 % Bit error rate, S2
switch LINK
    case "Uplink"                           
        % This is the TRACKING LOOP BANDWIDTH, NOT THE NOISE BANDWIDTH!
        CNR_MARGIN = 12;                            % dB
    case "Downlink"
        CNR_MARGIN = 10;
end

MODULATION = "BPSK";                                % S2
MODULATION_COEF = 1;                                % Class notes
MODULATION_INDEX = 72;                              % º, S2
ENCODING = "Turbo 1/6, 8920 frame length";          % S2
% They do use other encoding schemes like concatenated convolutional + RS
% but it's for lower bit rates
ENCODING_COEF = 6;

ANTENNA_T = 250; % deg, S2
ANTENNA_ACCURACY = 0.069;                                           % K, S2
switch ANTENNA % S2
    case "HGA"
        ANTENNA_D = 2.5;                            % deg, S2
        ANTENNA_BEAMWIDTH_REAL = 0.25;              % deg, S2
        ANTENNA_EFFICIENCY = 0.55; % Class notes for parabolic, even though this is a dual reflector        
        switch BAND
            case "X"
                switch LINK
                    case "Uplink"
                        ANTENNA_G_REAL = 43;
                    case "Downlink"
                        ANTENNA_G_REAL = 44.5;
                end
        end
    case "TLGA"
        ANTENNA_D = 0.26;
        ANTENNA_BEAMWIDTH_REAL = 10;
        switch LINK
            case "Uplink"
                ANTENNA_G_REAL = 5.5;
            case "Downlink"
                ANTENNA_G_REAL = 6.5;
        end
    case "MGA"
        ANTENNA_D = 0.135;
        ANTENNA_EFFICIENCY = 0.52;
        ANTENNA_BEAMWIDTH_REAL = 10;
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
        DSN_T = 21;                                 % K, S5
        switch LINK
            case "Uplink"
                DSN_G_REAL = 66.93;
                BANDWIDTH = 100;
            case "Downlink"
                DSN_G_REAL = 68.26;
                BANDWIDTH = 3;
        end
end

%% Loading the target coordinates https://ssd.jpl.nasa.gov/horizons/app.html#/
[Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno();

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
symbolRate = DATARATE * ENCODING_COEF / MODULATION_COEF;
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
        DSN_T = 49.72;
        N0receiver = 10*log10(BOLTZMANN*DSN_T)+30;
end
fprintf(strcat("Receiver noise spectral density: ",string(N0receiver)," dBm/Hz\n"));
fprintf(strcat("Received Pt/No: ",string(powerReceiver - N0receiver)," dB-Hz\n"));

% EBN0 - Error per Bit to Noise density
fprintf(strcat("Required Pt/No: ",string(10*log10(DATARATE))," dB-Hz\n"));
ebN0 = powerReceiver - N0receiver - 10*log10(DATARATE);
fprintf(strcat("Pt/No margin / Eb/No: ",string(ebN0)," dB\n"));

fprintf(strcat("Maximum bitrate for a minimum 0.7 dB (BER1e-5): ",string(10^((powerReceiver-0.7-N0receiver)/10)/1000),"kbps \n"));

% CARRIER MODULATION INDEX REDUCTION
powerModLoss = 20*log10(cosd(MODULATION_INDEX));
if (BAND=="Ka")
    powerModLoss = 0;
end
fprintf(strcat("Carrier Modulation Index Reduction: (telemetry carrier suppresion) ",string(powerModLoss)," dB\n"))

% CARRIER POWER
powerCarrier = powerReceiver + powerModLoss;
fprintf(strcat("Carrier power: ",string(powerCarrier)," dBm\n"))
% Pc/N0
powerCarrierNoiseRatio = powerCarrier-N0receiver;
fprintf(strcat("Pc/N0: ",string(powerCarrierNoiseRatio)," dB\n"))

% BANDWIDTH
fprintf(strcat("Bandwidth: ",string(BANDWIDTH)," Hz or ",string(10*log10(BANDWIDTH))," dB\n"))

% SIGNAL TO NOISE RATIO
snrCarrier = powerCarrierNoiseRatio - 10*log10(BANDWIDTH);
fprintf(strcat("CNR: ",string(snrCarrier)," dB\n"))

snrMargin = snrCarrier - CNR_MARGIN;
fprintf(strcat("CNR margin > ",string(CNR_MARGIN),": ",string(snrMargin)," dB\n"))

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

performance2.kbpsMaxforEbN0 = 10.^((PtN0-0.7)/10)/1000;

maxy = prctile(performance2.kbpsMaxforEbN0,80);
performance2.kbpsMaxforEbN0(performance2.kbpsMaxforEbN0 > maxy) = maxy;
performance3.Distance = distance.JunoToEarth;   
Pc = Prx + powerModLoss;
PcN0 = Pc - N0rx;
performance3.CNR = PcN0 - 10*log10(BANDWIDTH);
performance3.CNRthreshold = 3+ones(size(performance.Distance))*CNR_MARGIN;
cnrMargin = performance3.CNR - performance3.CNRthreshold;

plotOnEngine2D( ...
    time, ...
    performance, ...
    strcat("Losses for ",LINK,"-",ANTENNA), ...
    'Distance to Earth [AU]', ...
    'Loss [dB]', ...
    'Date [year]',true);

plotOnEngine2D( ...
    time, ...
    performance2, ...
    strcat("Max bitrate for ",LINK,"-",ANTENNA), ...
    'Distance to Earth [AU]', ...
    'Bitrate [kbps]', ...
    'Date [year]');

plotOnEngine2D( ...
    time, ...
    performance3, ...
    strcat("Carrier loop SNR for ",LINK,"-",ANTENNA), ...
    'Distance to Earth [AU]', ...
    'CNR [dB]', ...
    'Date [year]');