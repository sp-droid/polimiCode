function a = srp(r,rSun,TTsun,Rsun,cR,AoverM)

% Speed of light
c = 299792458;
% Stefan-Boltzmann constant in W/k^4/m^2
sigm = 5.67*1e-8;

% Surface solar flux
phi0 = sigm*TTsun^4;

% Solar flux near Earth ~ near the satellite because rSun >> r
phi = phi0*(Rsun/norm(rSun))^2;
solarPressure = phi/c;

dirSun = normalize(rSun,'norm');
a = -solarPressure*cR*AoverM*dirSun/1000;
end