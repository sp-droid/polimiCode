function a = drag(r, v, rho, wEarth, cD, AoverM)
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
%H0 = R(gas constant)*T/grav acceleration at that height/mean molecular
%weight of constituents
% vrel, cD and density rho are contentious parts of this approximation

vrel = (v - cross([0;0;wEarth],r));
%AoverM[m^2/kg], cD[-], rho[kg/m^3], vrel^2[km/s]^2 --> km^2/m/s^2 *1000-> km/s^2
a = (-0.5*AoverM*cD*rho*norm(vrel)*vrel)*1000;
end