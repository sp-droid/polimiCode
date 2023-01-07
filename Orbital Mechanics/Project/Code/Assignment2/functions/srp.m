function a = srp(r,rSun,TTsun,Rsun,Rbody,cR,AoverM)
% INPUTS:
% I had no time to write the descriptions, i'm leaving them like this for now. The scripts on Code/aasignment2 should be self explanatory, in particular the realcomparison.m. It's 100% updated
%
% OUTPUT:
% a[3x1] Perturbing acceleration
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------

% Speed of light
c = 299792458;
% Stefan-Boltzmann constant in W/k^4/m
sigm = 5.67*1e-8;

%% Solar flux
% Surface solar flux
phi0 = sigm*TTsun^4;

% Solar flux near Earth ~ near the satellite because rSun >> r
phi = phi0*(Rsun/norm(rSun))^2;
solarPressure = phi/c; % kg/m/s^2

%% Shadow function
% Apparent radius of the sun
a = asin(Rsun/norm(rSun-r));
% Apparent radius of the occulting body
b = asin(Rbody/norm(r));
% Apparent separation of centers
c = acos(-dot(r,rSun-r)/norm(r)/norm(rSun-r));
% Partial occultation
if (abs(a-b) < c) && (c < a+b)
    x = (c^2+a^2-b^2)/2/c;
    y = sqrt(a^2-x^2);
    
    A = a^2*acos(x/a) + b^2*acos((c-x)/b) - c*y;
    nu = 1-A/pi/a^2;
else
    % No occultation
    if (a+b < c)
        nu = 1;
    % Total occultation or partial but maximum
    else
        nu = 0;
    end
end

a = -nu*solarPressure*cR*AoverM*normalize(rSun,'norm')/1000;
end