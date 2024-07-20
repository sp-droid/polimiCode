%% Tracking Telemetry & Telecommand subsystem sizing
clc
clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.25 0.25 0.3 0.4]);
addpath('..\..\commonFunctions')

%% Loading the target coordinates https://ssd.jpl.nasa.gov/horizons/app.html#/
[Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno();

%% Parameters
ANTENNA = "HGA";
BAND = "X";
for j=1:2
    if j==1
        LINK = "Uplink";
    else
        LINK = "Downlink";
    end
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
    
    %% Constants
    AU = 149597870700;%m
    BOLTZMANN = 1.38e-23;%Ws/K
    
    %% Sizing
    
    % AMPLIFIER
    input_power = AMP_OUT_POWER / AMP_EFFICIENCY;
    
    % SYMBOL RATE
    symbolRate = DATARATE * ENCODING_COEF / MODULATION_COEF;
    
    wavelength = 3e8/FREQ_CARRIER;
    
    % JUNO ANTENNA
    antennaGain = 10*log10(pi^2 * ANTENNA_D^2 * ANTENNA_EFFICIENCY / wavelength^2);
    antennaBeamwidth = 65.3 * wavelength / ANTENNA_D; % https://en.wikipedia.org/wiki/Parabolic_antenna#Beamwidth

    % GROUND STATION ANTENNA
    dsnGain = 10*log10(pi^2 * DSN_D^2 * DSN_EFFICIENCY / wavelength^2);
    dsnBeamwidth = 65.3 * wavelength / DSN_D;
    
    % POINTING LOSSES
    % The pointing accuracy of the DSN antenna was calculated by taking the real pointing loss of -0.1 dB
    % The one from Juno was fine tuned but mostly kept the same as the recommended 0.1º
    
    lossesPointing = -12*(ANTENNA_ACCURACY/ANTENNA_BEAMWIDTH_REAL)^2;
    lossesPointing = lossesPointing + -12*(DSN_ACCURACY/dsnBeamwidth)^2;
    
    % ATMOSPHERIC LOSSES (ignoring rain)
    lossesAtm = -0.2; %This model is nice https://www.mathworks.com/help/phased/ref/gaspl.html, but we need to input height and temperature profiles...
    
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
    eirp = 10*log10(transmitterOutPower)+30 + transmitterGain + lossesCables;
    
    % FREE SPACE LOSSES - FROM THIS PART, THE DISTANCE MATTERS
    maxDistance = max(distance.JunoToEarth)*AU;
    lossesSpace = 20*log10(wavelength/4/pi/maxDistance);
    
    powerReceiver = eirp + receiverGain + lossesSpace + lossesPointing + lossesAtm;
    
    % SYSTEM NOISE DENSITY
    switch LINK
        case "Uplink"
            N0receiver = 10*log10(BOLTZMANN*ANTENNA_T)+30;
        case "Downlink"
            N0receiver = 10*log10(BOLTZMANN*DSN_T)+30;
    end
    
    % EBN0 - Error per Bit to Noise density
    ebN0 = powerReceiver - N0receiver - 10*log10(DATARATE);
    
    % CARRIER MODULATION INDEX REDUCTION
    powerModLoss = 20*log10(cosd(MODULATION_INDEX));
    
    % CARRIER POWER
    powerCarrier = powerReceiver + powerModLoss;
    % Pc/N0
    powerCarrierNoiseRatio = powerCarrier-N0receiver;
    
    % SIGNAL TO NOISE RATIO
    snrCarrier = powerCarrier - N0receiver - 10*log10(BANDWIDTH);
    
    snrMargin = snrCarrier - CNR_MARGIN;
    
    %% Plots
    totalLoss = 20*log10(wavelength/4/pi./distance.JunoToEarth/AU) + lossesPointing + lossesAtm + lossesCables;
    
    Prx = eirp + receiverGain + totalLoss - lossesCables;
    N0rx = N0receiver;
    PtN0 = Prx - N0rx;
    
    kbpsMaxforEbN0 = 10.^((PtN0-0.7)/10)/1000;
    maxy = prctile(kbpsMaxforEbN0,80);
    kbpsMaxforEbN0(kbpsMaxforEbN0 > maxy) = maxy;

    Pc = Prx + powerModLoss;
    PcN0 = Pc - N0rx;
    cnr = PcN0 - 10*log10(BANDWIDTH);
    CNRthreshold = 3+ones(size(totalLoss))*CNR_MARGIN;
    cnrMargin = cnr - CNRthreshold;

    if j==1
        performance.TotalLoss.Distance = distance.JunoToEarth;
        performance.kbpsMax.Distance = distance.JunoToEarth;
        performance.cnrMargin.Distance = distance.JunoToEarth;
    end
    performance.TotalLoss.(string(LINK)) = totalLoss;
    performance.kbpsMax.(string(LINK)) = kbpsMaxforEbN0;
    performance.cnrMargin.(string(LINK)) = cnrMargin;
end
plotOnEngine2D( ...
    time, ...
    performance.TotalLoss, ...
    strcat("Total losses for ",ANTENNA," antenna",""), ...
    '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
    '\boldmath$Loss\hspace{0.5em}[dB]$', ...
    '\boldmath$Date\hspace{0.5em}[year]$');

plotOnEngine2D( ...
    time, ...
    performance.kbpsMax, ...
    strcat("Max bitrate for ",ANTENNA," antenna",""), ...
    '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
    '\boldmath$Bitrate\hspace{0.5em}[kbps]$', ...
    '\boldmath$Date\hspace{0.5em}[year]$');

plotOnEngine2D( ...
    time, ...
    performance.cnrMargin, ...
    strcat("Carrier loop SNR margin for ",ANTENNA," antenna",""), ...
    '\boldmath$Distance\hspace{0.5em}to\hspace{0.5em}Earth\hspace{0.5em}[AU]$', ...
    '\boldmath$CNR\hspace{0.5em}margin\hspace{0.5em}[dB]$', ...
    '\boldmath$Date\hspace{0.5em}[year]$');