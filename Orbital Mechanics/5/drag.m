clear
close all

cD = 2.1;
Aoverm = 0.0095; %m^2/kg

% Assumptions:
% Ignoring lift and binormal aerodynamic forces
% Assuming no effect of attitude on cD
% Assuming constant cD, even though particles interact with the s/c
% in different ways depending on impact surface material, atmosphere
% chemistry, local air density & temperature...
% Wind direction parallel to s/c velocity
% Assuming air co-rotates with the Earth, though some very strong winds
% have been detected in the upper atmosphere (100s of m/s)
% Basic rho model
r = [ 7495.3;0.0000;0.0000 ];
v = [ 0.0000;0.2686;-7.3239 ];
%H0 = R(gas constant)*T/grav acceleration at that height/mean molecular
%weight of constituents
rho = rho0*exp(-(rnorm-Rearth)/H0);
% v = v - cross(wEarth, r)

a = -0.5*cD*Aoverm*rho*norm(v)*v