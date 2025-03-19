function att = pathAttenuationITUR676(elevAngle, freq, T, P, rhoW)
% Path attenuation model from ITU-R P.676-12 annex 2
% src: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.676-12-201908-S!!PDF-E.pdf
%
% PROTOTYPE
% att = attenuationITUR676( 1000, 10, 15, 101325, 7.5 )
%
% INPUTS:
% elevAngle[1]  - Elevation angle [deg]
% freq[1]	    - Frequency [GHz]
% T[1]	        - Ambient temperature [ºC]
% P[1]	        - Dry air pressure [Pa]
% rhoW[1]	    - Water vapor density [g/m^3]
%
% OUTPUTS:
% att[1]        - Attenuation [ dB ]
%
% NOTES:
% Good for frequencies [0,350] GHz and elevation angles [5,90] deg 
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------

att = attenuationITUR676( freq, T, P, rhoW, 1 );
attOx = att(1); attW = att(2);

% T in kelvin
T = T+273.15;
% P in hectopascals
P = P/100;
theta = 300/T;
% Water vapor partial pressure
e = rhoW*T/216.7;

rp = (P+e)/1013.25;

% Oxygen equivalent height
t1 = 5.1040/(1+0.066*rp^-2.3)*exp(-((freq-59.7)/(2.87+12.4*exp(-7.9*rp)))^2);
t2 = 0;
c = [0.1597 118.750334;0.1066 368.498246;0.1325 424.763020;0.1242 487.249273;0.0938 715.392902;0.1448 773.839490;0.1374 834.145546];
for i=1:length(c)
    ci = c(i,1); fi = c(i,2);
    t2 = t2 + ci*exp(2.12*rp)/((freq-fi)^2+0.025*exp(2.2*rp));
end
t3 = 0.0114*freq/(1+0.14*rp^-2.6)*(15.02*freq^2-1353*freq+5.333*1e4)/(freq^3-151.3*freq^2+9629*freq-6803);
A = 0.7832+0.00709*(T-273.15);
hOx = 6.1*A/(1+0.17*rp^-1.1)*(1+t1+t2+t3);

if freq<70 && hOx>10.7*rp^0.3
    hOx = 10.7*rp^0.3;
end

% Water equivalent height
sigmaW = 1.013/(1+exp(-8.6*(rp-0.57)));
ab = [22.235080 1.52 2.56;183.310087 7.62 10.2;325.152888 1.56 2.70;380.197353 4.15 5.70;439.150807 0.20 0.91;448.001085 1.63 2.46;...
    474.689092 0.76 2.22;488.490108 0.26 2.49;556.935985 7.81 10.0;620.70087 1.25 2.35;752.033113 16.2 20.0;916.171582 1.47 2.58;...
    970.315022 1.36 2.44;987.926764 1.60 1.86];
hW = 0;
for i=1:length(ab)
    fi = ab(i,1); ai = ab(i,2); bi = ab(i,3);
    hW = hW + ai*sigmaW/((freq-fi)^2+bi*sigmaW);
end
A = 1.9298-0.04166*(T-273.15)+0.0517*rhoW;
B = 1.1674-0.00622*(T-273.15)+0.0063*rhoW;
hW = A+B*hW;

Aox = attOx*hOx;
Aw = attW*hW;

att = (Aox+Aw)/sind(elevAngle);

end