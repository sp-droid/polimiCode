function rho = densityHigh(h, f107, Ap)
% INPUTS:
% h[1]      Height, in km. Only valid in high altitudes
% F107[1]   Solar 10.7 cm radio emissions, in solar flux units. [65,300]
% Ap[1]     Geomagnetic index. Usually [0,400]
% OUTPUT:
% rho[1]	Density 
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------

T = 900+2.5*(f107-70)+1.5*Ap;   % Temperature
m = 27-0.012*(h-200);           % Molecular mass
H = T/m;

rho = 6e-10*exp(-(h-175)/H);

end