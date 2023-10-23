function [kep] = car2kep (r,v, mu)
% 
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% DESCRIPTION:
%   Conversion from Cartesian coordinates to Keplerian elements
%
% INPUT
%   r       [3x1]   Position vector         [km]
%   v       [3x1]   Velocity vector         [km/s]
%   mu      [1x1]   Gravitational parameter [km^3/s^2]
%
% OUTPUT
%   kep     [1x6]   Vector of orbit keplerian elements. It contains:
%                       a       [1x1]   Semi-major axis         [km]
%                       e       [1x1]   Eccentricity            [-]
%                       i       [1x1]   Inclination             [rad]
%                       OM      [1x1]   RAAN                    [rad]
%                       om      [1x1]   Pericentre anomaly      [rad]
%                       th      [1x1]   True anomaly            [rad]
%
% VERSIONS:
%   2022-10-13
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

%% State and velocity vectors
r_norm = norm(r);
v_norm = norm(v);
%% Angular momentum h
h = cross(r,v);
h_norm = norm(h);
%% Inclination i
i = acos(h(3)/h_norm);
%% Eccentricity
e = 1/mu * ((v_norm^2 - mu/r_norm)*r - (dot(r,v))*v);
e_norm = norm(e);
%% Specific energy and semi-major axis
E = 0.5 * v_norm^2 - mu/r_norm;
a = - mu / (2*E);
%% Node line
k = [0;0;1];
if h(1) == 0 && h(2) == 0
    N = [1;0;0];
else
    N =  cross(k,h);
end
N_norm = norm(N);
%% Right ascension of ascending node
if N(2)>=0
    OM = acos(N(1)/N_norm);
else
    OM = 2*pi - acos(N(1)/N_norm);
end
if i == 0
    OM = 0;
end
%% Pericentre anomaly
if e(3)>=0
    if e_norm == 0
        om = 0;
    else
        om = acos(dot(N,e)/(N_norm*e_norm));
    end
else
    om = 2*pi - acos(dot(N,e)/(N_norm*e_norm));
end
%% Radial velocity
v_r = dot(r,v) / r_norm;
%% True anomaly
if v_r >= 0
    if e_norm == 0
        th = acos(dot(N,r)/(N_norm*r_norm));
    else
        th = acos(dot(e,r)/(e_norm*r_norm));
    end
else
    if e_norm == 0
        th = 2*pi - acos(dot(N,r)/(N_norm*r_norm));
    else
        th = 2*pi - acos(dot(e,r)/(e_norm*r_norm));
    end
end

kep = [a, e_norm, i, OM, om, th];










