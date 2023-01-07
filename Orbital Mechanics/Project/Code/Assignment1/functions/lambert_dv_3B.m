function dv = lambert_dv_3B(date,ID1,ID2,ID3,mu,opt_lambert_solv,h_atm)
% 
% This function computes the total cost needed to go from body 1 to body 3 with a flyby at body 2
% for given departure, fly-by and arrival dates.
% 
% INPUTS:
%   date    Vector containing in order departure, fly-by and arrival date
%           expressed in mjd2000
%   ID1     Natural number representing departure body 1
%   ID2     Natural number represeting fly-by body 2
%   ID3     Natural number represeting arrival body 3
%   mu      Planetary constant (Needed for bodies orbits and transfer
%           trajecotries)
%   h_atm   Atmosphere scale height of fly-by body
% 
% OUTPUT:
%   dv      Total cost of the given mission
% 
% VERSIONS:
%   2022-12-27: First version
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

orbitType1 = opt_lambert_solv.orbitType1; orbitType2 = opt_lambert_solv.orbitType2;

date1 = date(1);
date2 = date(2);
date3 = date(3);

t1 = date1 * (24*60*60);
t2 = date2 * (24*60*60);
t3 = date3 * (24*60*60);
tof1 = t2 - t1;
tof2 = t3 - t2;
% Orbits keplerian elements
if ID1 < 12
    [kep1,~] = uplanet(date1, ID1);
else
    [kep1,~,~] = ephNEO(date1,ID1);
end
if ID2 < 12
    [kep2,~] = uplanet(date2, ID2);
else
    [kep2,~,~] = ephNEO(date2, ID2);
end
if ID3 < 12
    [kep3,~] = uplanet(date3, ID3);
else
    [kep3,~,~] = ephNEO(date3, ID3);
end
% Orbits cartesian elements
[r1, v1] = kep2car(kep1, mu);
[r2, v2] = kep2car(kep2, mu);
[r3, v3] = kep2car(kep3, mu);
% Lambert solver
[~,~,~,ERR1,v1_t1,v2_t1,~,~] = lambertMR(r1, r2, tof1, mu, orbitType1, 0, 0, 2);
dv_1 = v1_t1' - v1;
[~,~,~,ERR2,v1_t2,v2_t2,~,~] = lambertMR(r2, r3, tof2, mu, orbitType2, 0, 0, 2);
dv_2 = v3 - v2_t2';
% Flyby
v_inf_m = v2_t1' - v2; v_inf_p = v1_t2' - v2; % s/c velocities relative to the planet
[r_p, ~, delta_v, ~, ~, ~] = power_gravity_assist(v_inf_m, v_inf_p, ID2, h_atm);
% Total cost
dv = NaN;
if all([~ERR1, ~ERR2, isfinite(r_p')])
    dv = abs(norm(dv_1)) + abs(norm(dv_2)) + abs(delta_v);
end