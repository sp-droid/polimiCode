function  [r_p, delta, delta_v, kep_i, kep_f, centres] = power_gravity_assist(v_i, v_f, ID, h_atm)
% 
% This function computes a powered gravity assist flyby around the specified planet: two hyperbolic arcs are
% combined by performing a tangential impulsive manoeuvre at their common pericentre.
% 
% INPUT:
%   v_i         [3x1]   Incoming velocity of first hyperbola at infinite
%   v_f         [3x1]   Outgoing velocity of second hyperbola at infinite
%   ID          [1x1]   Planet ID:
%                           1: Mercury          6: Saturn
%                           2: Venus            7: Uranus
%                           3: Earth            8: Neptune
%                           4: Mars             9: Pluto
%                           5: Jupiter          10: Moon
%   h_atm       [1x1]   Atmosphere height
% 
%   OUTPUT:
%   r_p         [3x1]   Radius of pericentre vector (NaN if the solution has not been found) 
%   delta       [1x1]   Turning angle
%   delta_v     [1x1]   Magnitude of impulsive manoeuvre given at common pericentre
%                       (NaN if the solution has not been found)
%   kep_i       [1x6]   Vector [a_i, e_i] with a_i and e_i respectively the
%                       semi-major axis and eccentricity of first hyperbola
%                       (NaN if the solution has not been found)
%   kep_f       [1x6]   Vector [a_f, e_f] with a_f and e_f respectively the
%                       semi-major axis and eccentricity of second hyperbola
%                       (NaN if the solution has not been found)
%   centres     [3x2]   Matrix containing the coordinates of the centres of the two 
%                       hyperbolae, respectively on the first and on the
%                       second column
%
% CONTRIBUTORS
%   Alessandro Michelazzi
%
% VERSIONS
%   2022-12-14: First version

% Data
if ID == 0
    % Sun
    R_p = astroConstants(3);
    mu = astroConstants(4);
elseif all([ID > 0, ID < 11])
        R_p = astroConstants(ID + 20);
        mu = astroConstants(ID + 10);
else
    error('Wrong ID. The function can only be used for flybys around Sun, planets and Moon')
end
if nargin < 4
    h_atm = 0;
end

R_p = R_p + h_atm; % Minimum feasible pericentre

v_i_norm = norm(v_i);
v_f_norm = norm(v_f);
v_i_unitvec = v_i / v_i_norm;

e_i_fun = @(r) 1 + (r * v_i_norm^2) / mu;
delta_i_fun = @(r) 2 * asin(1/e_i_fun(r));

e_f_fun = @(r) 1 + (r * v_f_norm^2) / mu;
delta_f_fun = @(r) 2 * asin(1/e_f_fun(r));

v_i_unit  = v_i / norm(v_i);
v_f_unit = v_f / norm(v_f);

% Turning angle
delta = acos(dot(v_i_unit, v_f_unit));

delta_fun = @(r) (delta_i_fun(r)/2 + delta_f_fun(r)/2 - delta);

r_p0 = R_p*1.001; % Initial guess
options = optimset('Display', 'off');
[r_p_norm,~, exitflag] = fzero(delta_fun, r_p0,options);
solutionFound = 1;

if r_p_norm <= R_p % Radius of pericentre less than the minimum possible
    r_p = NaN(3,1);
    solutionFound = 0;
end
if exitflag ~= 1 % Solution not found
    r_p = NaN(3,1);
    solutionFound = 0;
end

if solutionFound
    % Compute hyperbolae
    e_i = e_i_fun(r_p_norm);
    e_f = e_f_fun(r_p_norm);
    a_i = r_p_norm / (1-e_i);
    a_f = r_p_norm / (1-e_f);
    kep_i = [a_i, e_i];
    kep_f = [a_f, e_f];
    h_i = sqrt(mu * a_i * (1-e_i^2));
    h_f = sqrt(mu * a_f * (1-e_f^2));
    v_p_i = mu / h_i * (1 + e_i);
    v_p_f = mu / h_f * (1 + e_f);
    delta_v = v_p_f - v_p_i;
    
    % Vector of radius of pericentre
    beta = pi/2 - delta_i_fun(r_p_norm)/2;
    u_unitvec = cross(v_i, v_f) / norm(cross(v_i, v_f));
    r_p_unitvec = v_i_unitvec*cos(-beta) + cross(u_unitvec,v_i_unitvec)*sin(-beta) ...
                + u_unitvec*dot(u_unitvec,v_i_unitvec)*(1-cos(-beta));
    r_p = r_p_norm * r_p_unitvec;
    
    % Centres of hyperbolae
    centre1_norm = r_p_norm + abs(a_i); 
    centre1 = centre1_norm * r_p_unitvec;
    centre2_norm = r_p_norm + abs(a_f); centre2 = centre2_norm * r_p_unitvec;
    centres = [centre1, centre2];
else
    delta_v = NaN(1); kep_i = NaN(1,2); kep_f = NaN(1,2); centres = NaN(3,2);
end





