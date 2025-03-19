function att = attenuationITUR676(freq, T, P, rhoW, outputType)
% Attenuation model from ITU-R P.676-12 annex 1
% src: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.676-12-201908-S!!PDF-E.pdf
%
% PROTOTYPE
% att = attenuationITUR676( 1000, 10, 15, 101325, 7.5 )
%
% INPUTS:
% freq[1]	    - Frequency [GHz]
% T[1]	        - Ambient temperature [ºC]
% P[1]	        - Dry air pressure [Pa]
% rhoW[1]	    - Water vapor density [km]
%
% OPTIONAL INPUTS:
% outputType[1] - 0 to return attenuation, 1 to return [oxygen,vapor] atts
%
% OUTPUTS:
% att[1]        - Specific attenuation [ dB/km ]
%
% NOTES:
% Good for frequencies [0,1000] GHz
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------

if nargin < 5
    outputType = 0;
end

% Load coefficients
load("v12OxygenLines.mat");
load("v12VaporLines.mat");

% T in kelvin
T = T+273.15;
% P in hectopascals
P = P/100;
theta = 300/T;
% Water vapor partial pressure
e = rhoW*T/216.7;

%%

% Imaginary part of freq-dependent complex refractivity
% Oxygen
NppOx = 0; 
for i=1:length(v12linesoxygen)
    fi = v12linesoxygen(i,1); ai = v12linesoxygen(i,2:7);

    delta = (ai(5) + ai(6)*theta)*1e-4*(P+e)*theta^0.8;
    deltaFreq = ai(3)*1e-4*(P*theta^(0.8-ai(4)) + 1.1*e*theta);
    deltaFreq = sqrt(deltaFreq^2+2.25*1e-6);
    Fi = freq/fi*((deltaFreq-delta*(fi-freq))/((fi-freq)^2+deltaFreq^2) + (deltaFreq-delta*(fi+freq))/((fi+freq)^2+deltaFreq^2));
    Si = ai(1)*1e-7 * P * theta^3 * exp(ai(2)*(1-theta));

    NppOx = NppOx + Si*Fi;
end
% Debye spectrum width parameter
d = 5.6*1e-4*(P+e)*theta^0.8;
% Dry continuum due to pressure-induced nitrogen absorption and Debye
% spectrum
NppD = freq*P*theta^2*(6.14*1e-5/d/(1+(freq/d)^2) + 1.4*1e-12*P*theta^1.5/(1+1.9*1e-5*freq^1.5));
NppOx = NppOx + NppD;

% Water vapor
NppW = 0;
for i=1:length(v12lineswatervapour)
    fi = v12lineswatervapour(i,1); bi = v12lineswatervapour(i,2:7);

    delta = 0;
    deltaFreq = bi(3)*1e-4*(P*theta^bi(4) + bi(5)*e*theta^bi(6));
    deltaFreq = 0.535*deltaFreq + sqrt(0.217*deltaFreq^2+1/theta*2.1316*1e-12*fi^2);
    Fi = freq/fi*((deltaFreq-delta*(fi-freq))/((fi-freq)^2+deltaFreq^2) + (deltaFreq-delta*(fi+freq))/((fi+freq)^2+deltaFreq^2));
    Si = bi(1)*0.1 * e * theta^3.5 * exp(bi(2)*(1-theta));
    NppW = NppW + Si*Fi;
end

% Gas attenuation = oxygen + water vapour attenuation
if outputType==0
    att = 0.182*freq*(NppOx+NppW);
else
    att = 0.182*freq*[NppOx, NppW];
end

end