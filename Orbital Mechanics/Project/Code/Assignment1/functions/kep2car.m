function [r,v] = kep2car(kep, mu)
% 
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE
%   [r, v] = kep2car(kep, mu)
%
% DESCRIPTION
%   Conversion from Keplerian elements to Cartesian coordinates Angles in radians
%
% INPUT
%   kep     [1x6]   Vector of orbit keplerian elements. It contains:
%                       a       [1x1]   Semi-major axis         [km]
%                       e       [1x1]   Eccentricity            [-]
%                       i       [1x1]   Inclination             [rad]
%                       OM      [1x1]   RAAN                    [rad]
%                       om      [1x1]   Pericentre anomaly      [rad]
%                       th      [1x1]   True anomaly            [rad]
%   mu      [1x1]   Gravitational parameter [km^3/s^2]
%
% OUTPUT
%   r       [3x1]   Position vector         [km]
%   v       [3x1]   Velocity vector         [km/s]
% 
% VERSIONS:
%   2022-10-12
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

a = kep(1); e = kep(2); i = kep(3);
OM = kep(4); om = kep(5); th = kep(6);

if th < 0
    th = th + 2*pi;
end

%% p
p = a * (1 - e^2);
%% r
r_norm = p / (1 + e*cos(th));
%% State vector
r_pf = r_norm .* [cos(th); sin(th); 0];
v_pf = sqrt(mu/p) .* [-sin(th); e+cos(th); 0];
%% Rotational matrices
R3_OM = [cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];
R1_i = [1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];
R3_om = [cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];

T = R3_OM' * R1_i' * R3_om';

%% State and velocity vectori in Inertial Earth Centred reference frame
r_ECI = T * r_pf;
v_ECI = T * v_pf;

r = r_ECI;
v = v_ECI;








